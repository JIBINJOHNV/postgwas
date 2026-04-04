import polars as pl
import pandas as pd
import pyarrow,sys
import gzip, lzma, zipfile, csv, subprocess, os,json,shlex,argparse
from typing import Tuple
from io import BytesIO,TextIOWrapper
from pathlib import Path
import yaml
from importlib import resources
from postgwas.harmonisation.chr_pos_process import fix_chr_pos_column
import os
import json

# ----------------------------------------------------------------------
# 1. Helper: correct opener
# ----------------------------------------------------------------------
def _opener(path: str):
    p = Path(path)
    if p.suffix.lower() == ".gz":      return gzip.open
    if p.suffix.lower() in {".xz", ".lzma"}: return lzma.open
    if p.suffix.lower() == ".zip":     return None 
    return open


# ----------------------------------------------------------------------
# 2. Delimiter detection – look at the *header* line
# ----------------------------------------------------------------------
def detect_delimiter(path: str) -> str:
    opener = _opener(path)
    if opener is None: # Zip
        try:
            with zipfile.ZipFile(path) as zf:
                csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
                if not csvs: return "\t"
                txt = zf.read(csvs[0]).decode(errors="ignore")
                header = next((l for l in txt.splitlines() if not l.startswith("##")), "")
        except:
            return "\t"
    else:
        try:
            with opener(path, "rt", errors="ignore") as f:
                header = next((l for l in f if not l.startswith("##")), "")
        except:
            return "\t"
    
    if not header.strip(): return "\t"
    try:
        return csv.Sniffer().sniff(header, delimiters="\t,; ").delimiter
    except csv.Error:
        return "\t"


# ----------------------------------------------------------------------
# 3. Detect ## metadata
# ----------------------------------------------------------------------
def has_double_hash(path: str) -> bool:
    opener = _opener(path)
    if opener is None: # Zip
        try:
            with zipfile.ZipFile(path) as zf:
                csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
                if not csvs: return False
                txt = zf.read(csvs[0]).decode(errors="ignore")
                return any(l.startswith("##") for l in txt.splitlines()[:50])
        except:
            return False
    try:
        with opener(path, "rt", errors="ignore") as f:
            for _ in range(50):
                line = f.readline()
                if not line: break
                if line.startswith("##"): return True
        return False
    except:
        return False


# ----------------------------------------------------------------------
# 4. Count data lines – **safe** shell command (list form)
# ----------------------------------------------------------------------
def count_data_lines(path: str, skip_hash: bool) -> int:
    if Path(path).suffix.lower() == ".zip":
        try:
            with zipfile.ZipFile(path) as zf:
                csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
                if not csvs: return 0
                txt = zf.read(csvs[0]).decode(errors="ignore")
                lines = txt.splitlines()
                if skip_hash:
                    lines = [l for l in lines if not l.startswith("##")]
                return len(lines)
        except Exception as e:
            print(f"Zip count failed: {e}")
            return 0
    # macOS-safe: use "< path" for compressed files
    if path.endswith(".gz"):
        decomp = "zcat <"
    elif path.endswith((".xz", ".lzma")):
        decomp = "xzcat <"
    else:
        decomp = "cat"
    cmd_str = f"{decomp} {shlex.quote(path)}"
    if skip_hash:
        cmd_str += " | grep -v '^##'"
    cmd_str += " | wc -l"
    try:
        result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True, check=True)
        return int(result.stdout.strip())
    except Exception as e:
        print(f"Shell count failed ({e}) — using Polars count only")
        return 0



# ----------------------------------------------------------------------
# 5. Helper: Clean Stream + Count
# ----------------------------------------------------------------------
def vcf_to_polars_stream(path: str) -> Tuple[BytesIO, int]:
    """
    1. Removes ## metadata lines.
    2. Counts valid lines on the fly.
    """
    buf = BytesIO()
    row_count = 0
    opener = _opener(path) 
    if opener is None and not path.endswith('.zip'): opener = open
    try:
        with opener(path, "rt", errors="ignore") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                buf.write(line.encode("utf-8"))
                row_count += 1
    except Exception as e:
        print(f"   ❌ Error streaming file: {e}")
        return BytesIO(), 0
    buf.seek(0)
    data_count = max(0, row_count - 1)
    return buf, data_count


# ======================================================================
# HELPER 1: DataFrame Cleanup (Headers, Types, Chromosome)
# ======================================================================
def clean_dataframe(df: pl.DataFrame, chr_col: str) -> pl.DataFrame:
    """
    Applied immediately after ANY read (Attempt 1 or Attempt 2).
    Trims whitespace, fixes types, and standardizes chromosome column.
    """
    if df.height == 0: return df
    # 1. Clean Headers
    df = df.rename({c: c.strip() for c in df.columns})
    # 2. Clean String Values
    for col in df.columns:
        if df[col].dtype in (pl.Utf8, pl.String):
            df = df.with_columns(pl.col(col).str.strip_chars())
    # 3. Type Conversion
    valid_cols = []
    for col in df.columns:
        if col == chr_col:
            valid_cols.append(pl.col(col)) # Skip casting for CHR
            continue
        numeric_series = df[col].cast(pl.Float64, strict=False)
        if numeric_series.null_count() > df[col].null_count():
            valid_cols.append(pl.col(col)) # Keep original if casting fails
        else:
            valid_cols.append(numeric_series)
    df = df.with_columns(valid_cols)
    # 4. Chromosome Fix (12.0 -> 12)
    if chr_col and chr_col in df.columns:
        df = df.with_columns(
            pl.coalesce(
                pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String),
                pl.col(chr_col)
            ).alias(chr_col)
        )
    return df

# ======================================================================
# HELPER 2: Validation (Is the file broken?)
# ======================================================================
def is_bad_parsing(df: pl.DataFrame, expected_cols: int) -> bool:
    if df.width < 2 and expected_cols > 1: return True # Swallowed
    if df.width > (expected_cols + 5): return True     # Ghost columns
    if df.height > 0:
        null_ratio = df.null_count().sum_horizontal().sum() / (df.height * df.width)
        if null_ratio > 0.5: return True # Null flood
    return False

# ======================================================================
# HELPER 3: Repair (Linux Shell Fix)
# ======================================================================
def normalize_whitespace_file(input_path: str, output_dir: str) -> str:
    base_name = Path(input_path).name
    while Path(base_name).suffix: base_name = Path(base_name).stem
    clean_path = os.path.join(output_dir, f"{base_name}_cleaned.tsv")
    if input_path.lower().endswith(".gz"): cat_cmd = "gzip -dc"
    elif input_path.lower().endswith(".zip"): cat_cmd = "unzip -p"
    elif input_path.lower().endswith((".xz", ".lzma")): cat_cmd = "xz -dc"
    else: cat_cmd = "cat"
    # Fix: zcat -> remove ## -> squeeze spaces -> tab
    cmd = (
        f"{cat_cmd} {shlex.quote(input_path)} | "
        f"grep -v '^##' | "
        f"tr -s '[:blank:]' '\\t' > {shlex.quote(clean_path)}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
        if os.path.exists(clean_path) and os.path.getsize(clean_path) > 0:
            return clean_path
    except: pass
    return ""

# ======================================================================
# HELPER 4: Header Counter
# ======================================================================
def get_expected_cols(path: str) -> int:
    try:
        # (Same opener logic as before) ...
        opener = gzip.open if path.lower().endswith((".gz", ".bgz")) else open
        if path.lower().endswith(".zip"): return 0
        
        with opener(path, "rt", errors="ignore") as f:
            for line in f:
                if not line.startswith("##") and line.strip():
                    # FIX: Replace commas (and semicolons) with spaces first
                    clean_line = line.replace(",", " ").replace(";", " ").replace("\t", " ")
                    return len(clean_line.split())
    except: 
        return 0
    return 0


def clean_and_optimize_dataframe(df: pl.DataFrame) -> pl.DataFrame:
    """
    1. Trims whitespace from ALL string columns.
    2. Replaces empty strings "" with null (optional, but good practice).
    3. Attempts to convert string columns that look like numbers into actual Integers/Floats.
    """
    return (
        df.lazy()
        # Step 1: Trim whitespace from all string columns
        .with_columns(
            pl.col(pl.String).str.strip_chars()
        )
        # Step 2 (Optional): Convert empty strings to null
        .with_columns(
            pl.col(pl.String).replace("", None)
        )
        # Step 3: Collect immediately to perform type inference
        .collect()
        # Step 4: Downcast types (e.g., String "123" -> Int64)
        .select(pl.all().shrink_dtype())
    )

# ======================================================================
# MAIN FUNCTION
# ======================================================================
def read_sumstats(sumstat_file: str, output_dir: str, sample_column_dict: dict) -> Tuple[pl.DataFrame, int, int]:

    chr_col = sample_column_dict.get("chr_col")
    chr_pos_col = sample_column_dict.get("chr_pos_col")
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------
    # Pre-flight Checks
    # ------------------------------
    delim = detect_delimiter(sumstat_file)
    skip_hash = has_double_hash(sumstat_file)
    expected_cols = get_expected_cols(sumstat_file)

    df = pl.DataFrame()
    shell_cnt = 0

    # ------------------------------
    # ATTEMPT 1: Optimistic Read
    # ------------------------------
    try:
        if Path(sumstat_file).suffix.lower() == ".zip":
            with zipfile.ZipFile(sumstat_file) as zf:
                csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
                if not csvs:
                    raise ValueError("No CSV/TSV in .zip")

                txt = zf.read(csvs[0]).decode(errors="ignore")
                lines = txt.splitlines()

                if skip_hash:
                    lines = [l for l in lines if not l.startswith("##")]

                shell_cnt = max(0, len(lines) - 1)

                df = pl.read_csv(
                    BytesIO("\n".join(lines).encode("utf-8")),
                    separator=delim,
                    has_header=True,
                    quote_char=None,
                    infer_schema_length=0,
                    ignore_errors=True,
                    null_values=["NA", "na", ".", ""],
                    truncate_ragged_lines=True,
                )

        else:
            if skip_hash:
                source, shell_cnt = vcf_to_polars_stream(sumstat_file)
                df = pl.read_csv(
                    source,
                    separator=delim,
                    has_header=True,
                    quote_char=None,
                    infer_schema_length=0,
                    ignore_errors=True,
                    null_values=["NA", "na", ".", ""],
                    truncate_ragged_lines=True,
                )
            else:
                shell_cnt = count_data_lines(sumstat_file, skip_hash=False)
                if shell_cnt >= 1:
                    shell_cnt -= 1

                df = pl.read_csv(
                    sumstat_file,
                    separator=delim,
                    has_header=True,
                    quote_char=None,
                    infer_schema_length=0,
                    ignore_errors=True,
                    null_values=["NA", "na", ".", ""],
                    truncate_ragged_lines=True,
                )
        df = clean_and_optimize_dataframe(df)

        df, sample_column_dict = fix_chr_pos_column(
            chromosome="All_Chrs",
            df=df,
            sample_column_dict=sample_column_dict,
            drop_mt=True
        )

        df = clean_dataframe(df, chr_col)

    except Exception as e:
        print(f"\t\t\t ⚠️ Attempt 1 Error: {e}")
        df = pl.DataFrame()

    # ------------------------------
    # ATTEMPT 2: Repair
    # ------------------------------
    if is_bad_parsing(df, expected_cols):

        print(f"\t\t\t⚠️  Standard parsing failed (Rows: {df.height}, Cols: {df.width}). Switching to Repair Mode...")

        clean_file = normalize_whitespace_file(sumstat_file, output_dir)

        if clean_file:
            print(f"\t\t\t🔄 Reading repaired file: {clean_file}")

            try:
                df = pl.read_csv(
                    clean_file,
                    separator="\t",
                    has_header=True,
                    null_values=["NA", "na", ".", ""],
                    truncate_ragged_lines=True
                )

                df, sample_column_dict = fix_chr_pos_column(
                    chromosome="All_Chrs",
                    df=df,
                    sample_column_dict=sample_column_dict,
                    drop_mt=True
                )

                df = clean_dataframe(df, chr_col)

                shell_cnt = df.height

                print(f"\t\t\t✅ Recovery Successful: {df.height} rows loaded.")

            except Exception as e:
                print(f"❌ Failed to read repaired file: {e}")
        else:
            print("❌ Repair failed (could not create clean file).")

    # ------------------------------
    # FINAL COLUMN FIXES
    # ------------------------------
    chr_col = sample_column_dict.get("chr_col")
    pos_col = sample_column_dict.get("pos_col")

    # 🔹 Fix Chromosome (chr-safe, X/Y safe)
    if chr_col in df.columns:
        df = df.with_columns(
            pl.col(chr_col)
            .cast(pl.Utf8, strict=False)
            .fill_null("")
            .str.strip_chars()
            .str.replace(r"(?i)^chr", "", literal=False)
            .str.to_uppercase()
            .alias(chr_col)
        )

        df = df.with_columns(
            pl.when(pl.col(chr_col).is_in(["X", "23", "23.0"])).then(pl.lit("X"))
            .when(pl.col(chr_col).is_in(["Y", "24", "24.0"])).then(pl.lit("Y"))
            .when(pl.col(chr_col).is_in(["MT", "M", "25", "25.0"])).then(pl.lit("MT"))
            .otherwise(pl.col(chr_col))
            .alias(chr_col)
        )

        # 🔴 Validate chromosomes
        valid_chr = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

        invalid_chr_df = df.filter(~pl.col(chr_col).is_in(valid_chr))

        if invalid_chr_df.height > 0:
            print("❌ Invalid chromosome values detected:")
            print(invalid_chr_df.select(chr_col).unique())

            df = df.filter(pl.col(chr_col).is_in(valid_chr))

        # 🔴 Optional (recommended): autosomes only
        # AUTOSOMES_ONLY = True
        # if AUTOSOMES_ONLY:
        #     df = df.filter(pl.col(chr_col).is_in([str(i) for i in range(1, 23)]))

    if pos_col in df.columns:
        df = df.with_columns(
            pl.col(pos_col)
            .cast(pl.Float64, strict=False)
            .cast(pl.Int64, strict=False)
        )

    # ------------------------------
    # FINAL QC
    # ------------------------------
    if df.height < 5:
        print(f"⚠️  WARNING: File is empty or has too few rows ({df.height}). Exiting.")
        return df, shell_cnt, df.height

    if abs(shell_cnt - df.height) > 1:
        print(f"⚠️ Warning: Mismatch between line-count ({shell_cnt}) and Polars rows ({df.height}).")
    
    indent = "\t" * 3
    print(f"{indent} First 5 rows of sumstats {sumstat_file}:")
    preview_str = str(df.head())
    preview_str = "\n".join(indent + line for line in preview_str.splitlines())
    print(preview_str)
    print(f"{indent}")
    print(f"{indent} Columns of sumstats {sumstat_file}:")
    cols = ", ".join(df.columns)
    print(f"{indent}[COLUMNS] ({len(df.columns)} total):\n{indent}{cols}") 
    return df, shell_cnt, df.height, sample_column_dict



def check_file_truncation(file_path):
    # -------------------------
    # 1. HARD FAIL
    # -------------------------
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"❌ FILE NOT FOUND: {file_path}")

    if os.path.getsize(file_path) == 0:
        raise ValueError(f"❌ EMPTY FILE: {file_path}")

    warnings = []

    # -------------------------
    # 2. GZIP integrity ONLY
    # -------------------------
    if file_path.endswith(".gz"):
        if subprocess.run(["gzip", "-t", file_path]).returncode != 0:
            warnings.append("GZIP corrupted or truncated")

    # -------------------------
    # 3. VERY LIGHT tail check
    # -------------------------
    try:
        cmd = f"zcat {shlex.quote(file_path)} | tail -n 1" if file_path.endswith(".gz") else f"tail -n 1 {shlex.quote(file_path)}"
        last_line = subprocess.check_output(cmd, shell=True).decode().strip()

        # ONLY strong signal → abrupt ending
        if last_line.endswith("."):
            warnings.append("Last line ends abruptly → possible truncation")

    except:
        pass  # ignore silently (minimal philosophy)

    # -------------------------
    # 4. PRINT only if needed
    # -------------------------
    if warnings:
        print("\n⚠️ Possible file truncation detected:")
        for w in warnings:
            print(f"   - {w}")

    return warnings

# def read_sumstats(sumstat_file: str, output_dir: str, sample_column_dict: dict) -> Tuple[pl.DataFrame, int, int]:
#     chr_col = sample_column_dict.get("chr_col")
#     chr_pos_col = sample_column_dict.get("chr_pos_col")
#     os.makedirs(output_dir, exist_ok=True)
#     # Pre-flight Checks
#     delim = detect_delimiter(sumstat_file)
#     skip_hash = has_double_hash(sumstat_file)
#     expected_cols = get_expected_cols(sumstat_file)
#     df = pl.DataFrame()
#     shell_cnt = 0
#     # ------------------------------------------------------------------
#     # ATTEMPT 1: Optimistic Read (YOUR ORIGINAL LOGIC)
#     # ------------------------------------------------------------------
#     try:
#         if Path(sumstat_file).suffix.lower() == ".zip":
#             # ZIP Handling
#             with zipfile.ZipFile(sumstat_file) as zf:
#                 csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
#                 if not csvs: raise ValueError("No CSV/TSV in .zip")
#                 txt = zf.read(csvs[0]).decode(errors="ignore")
#                 lines = txt.splitlines()
#                 if skip_hash:
#                     lines = [l for l in lines if not l.startswith("##")]
#                 shell_cnt = max(0, len(lines) - 1)
#                 clean_txt = "\n".join(lines)
#                 df = pl.read_csv(
#                     BytesIO(clean_txt.encode("utf-8")),
#                     separator=delim,
#                     has_header=True,
#                     quote_char=None,
#                     infer_schema_length=0,
#                     ignore_errors=True,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True,
#                 )
#         else:
#             # GZ / XZ / Text Handling
#             if skip_hash:
#                 source, shell_cnt = vcf_to_polars_stream(sumstat_file)
#                 df = pl.read_csv(
#                     source,
#                     separator=delim,
#                     has_header=True,
#                     quote_char=None,
#                     infer_schema_length=0,
#                     ignore_errors=True,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True,
#                 )
#             else:
#                 shell_cnt = count_data_lines(sumstat_file, skip_hash=False)
#                 if shell_cnt >= 1: shell_cnt -= 1
#                 df = pl.read_csv(
#                     sumstat_file,
#                     separator=delim,
#                     has_header=True,
#                     quote_char=None,
#                     infer_schema_length=0,
#                     ignore_errors=True,
#                     null_values=["NA", "na", ".", ""],
#                     truncate_ragged_lines=True,
#                 )
#         df=clean_and_optimize_dataframe(df)
#         df,sample_column_dict=fix_chr_pos_column( chromosome="All_Chrs",
#             df=df,sample_column_dict=sample_column_dict,drop_mt=True)
#         df = clean_dataframe(df, chr_col)
#     except Exception as e:
#         print(f"\t\t\t ⚠️ Attempt 1 Error: {e}")
#         df = pl.DataFrame()
#     # ------------------------------------------------------------------
#     # ATTEMPT 2: Validation & Repair (If Attempt 1 was messy)
#     # ------------------------------------------------------------------
#     if is_bad_parsing(df, expected_cols):
#         print(f"\t\t\t⚠️  Standard parsing failed (Rows: {df.height}, Cols: {df.width}). Switching to Repair Mode...")
#         # 1. Run Repair Tool
#         clean_file = normalize_whitespace_file(sumstat_file, output_dir)
#         if clean_file:
#             print(f"\t\t\t🔄 Reading repaired file: {clean_file}")
#             try:
#                 # 2. Re-read (Force Tab Separator)
#                 df = pl.read_csv(
#                     clean_file, 
#                     separator="\t", 
#                     has_header=True, 
#                     null_values=["NA", "na", ".", ""], 
#                     truncate_ragged_lines=True
#                 )
#                 # 3. Cleanup Again (Clean Headers/Types on new DF)
#                 df,sample_column_dict=fix_chr_pos_column( chromosome="All_Chrs",
#                     df=df,sample_column_dict=sample_column_dict,drop_mt=True)
#                 df = clean_dataframe(df, chr_col)
#                 # 4. Update Counts
#                 shell_cnt = df.height
#                 print(f"\t\t\t✅ Recovery Successful: {df.height} rows loaded.")
#             except Exception as e:
#                 print(f"❌ Failed to read repaired file: {e}")
#         else:
#             print("❌ Repair failed (could not create clean file).")
#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     # Filter to ensure we only touch columns that actually exist in this dataframe
#     cols_to_fix = []
#     # # 1. Fix Chromosome: Float(12.0) -> Int(12) -> Str("12")
#     if chr_col in df.columns:
#         cols_to_fix.append(
#             pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String)
#         )
#     # # 2. Fix Position: Float(1.5e7) -> Int(15000000)
#     if pos_col in df.columns:
#         cols_to_fix.append(
#             pl.col(pos_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False)
#         )
#     # # 3. Apply all fixes at once
#     if cols_to_fix:
#         df = df.with_columns(cols_to_fix)
#     # ------------------------------------------------------------------
#     # FINAL QC & RETURN
#     # ------------------------------------------------------------------
#     # Stop if empty (checking after all attempts)
#     if df.height < 5:
#         print(f"⚠️  WARNING: File is empty or has too few rows ({df.height}). Exiting.")
#         return df, shell_cnt, df.height
#     # QC Warning
#     if abs(shell_cnt - df.height) > 1:
#         print(f"⚠️ Warning: Mismatch between line-count ({shell_cnt}) and Polars rows ({df.height}).")
#     return df, shell_cnt, df.height,sample_column_dict




def find_resource_file_path(input_file, resource_folder, grch_version, chromosome):
    """
    Resolve user EAF resource path.

    Logic:
      1. If input_file == "NA" -> return "NA"
      2. If input_file is a full existing path -> return it
      3. Otherwise build:
           {resource_folder}/{grch_version}/external_af/vcf_files/
           {grch_version}_{input_file}_freq_chr{chromosome}.vcf.gz
         If exists -> return it
      4. Else -> return "NA"
    """

    # 1) User explicitly disabled EAF
    if input_file == "NA":
        return "NA"

    # 2) User gave a full path
    if os.path.exists(input_file):
        return input_file
    
    if os.path.exists(f"{input_file}_chr{chromosome}.tsv.gz"):
        return f"{input_file}_chr{chromosome}.tsv.gz"
    
    # 3) Construct the standard resource path (EXACT pattern you gave)
    constructed_path = (
        f"{resource_folder}/{grch_version}/external_af/vcf_files/"
        f"{grch_version}_{input_file}_freq_chr{chromosome}.vcf.gz"
    )

    if os.path.exists(constructed_path):
        return constructed_path

    # 4) Nothing found
    print(f"⚠️ Warning: No valid EAF file found for '{input_file}'.")
    print(f"   Tried full path      : {input_file}")
    print(f"   Tried constructed path: {input_file}_chr{chromosome}.tsv.gz")
    return "NA"



def read_config(csv_path):
    """
    Read a configuration file (CSV or TSV) and return a list of dictionaries.
    Automatically detects delimiter and cleans whitespace.

    Parameters
    ----------
    csv_path : str
        Path to configuration file.

    Returns
    -------
    list[dict]
        List of configuration dictionaries.
    """

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"❌ The file '{csv_path}' was not found.")

    # -------------------------------
    # Auto-detect separator
    # -------------------------------
    df = pd.read_csv(csv_path, sep=None, engine="python")

    # Replace NA values
    df = df.fillna("NA")

    # -------------------------------
    # Clean column names
    # -------------------------------
    df.columns = df.columns.str.strip()

    # -------------------------------
    # Clean string cells
    # -------------------------------
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # -------------------------------
    # Convert to list of dictionaries
    # -------------------------------
    cfg_list = df.to_dict(orient="records")

    # Final safety strip
    cfg_list = [
        {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in cfg.items()}
        for cfg in cfg_list
    ]

    return cfg_list

def load_default_config(filename: str):
    """
    Load a YAML config file.

    New behaviour:
        • If `filename` is a valid filesystem path → load directly.
        • Else → try loading from postgwas.config package resource.
    """

    # 1️⃣ Try direct filesystem path
    if os.path.exists(filename) and os.path.isfile(filename):
        try:
            with open(filename, "r") as f:
                return yaml.safe_load(f)
        except Exception as e:
            raise RuntimeError(
                f"❌ ERROR: Failed to load config from path:\n  {filename}\nReason: {e}"
            )

    # 2️⃣ Fallback to package resource
    try:
        with resources.open_text("postgwas.config", filename) as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"❌ Config file '{filename}' not found.\n"
            f"Tried filesystem path AND postgwas.config resource."
        )
    except Exception as e:
        raise RuntimeError(
            f"❌ ERROR: Failed to load YAML from package resource '{filename}'.\nReason: {e}"
        )
