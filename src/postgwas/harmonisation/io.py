import polars as pl
import pandas as pd
import pyarrow
import gzip, lzma, zipfile, csv, subprocess, os,json
from typing import Tuple
from io import BytesIO
from pathlib import Path
import shlex
import argparse
import sys
import yaml
from importlib import resources

# ----------------------------------------------------------------------
# 1. Helper: correct opener
# ----------------------------------------------------------------------
def _opener(path: str):
    p = Path(path)
    if p.suffix.lower() == ".gz":      return gzip.open
    if p.suffix.lower() in {".xz", ".lzma"}: return lzma.open
    if p.suffix.lower() == ".zip":     return None          # handled separately
    return open


# ----------------------------------------------------------------------
# 2. Delimiter detection – look at the *header* line
# ----------------------------------------------------------------------
def detect_delimiter(path: str) -> str:
    opener = _opener(path)
    if opener is None:
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
    if not header.strip():
        return "\t"
    try:
        # csv.Sniffer works best on a *single* line
        return csv.Sniffer().sniff(header, delimiters="\t,; ").delimiter
    except csv.Error:
        return "\t"


# ----------------------------------------------------------------------
# 3. Detect ## metadata
# ----------------------------------------------------------------------
def has_double_hash(path: str) -> bool:
    opener = _opener(path)
    if opener is None:                     # .zip
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
# 5. Main function
# ----------------------------------------------------------------------
def read_sumstats(sumstat_file: str,output_dir: str) -> Tuple[pl.DataFrame, int, int]:
    os.makedirs(output_dir, exist_ok=True)
    #print(f"Reading: {sumstat_file}")
    # ---- delimiter ----------------------------------------------------
    delim = detect_delimiter(sumstat_file)
    #print(f"Delimiter: '{delim}'")
    # ---- ## metadata --------------------------------------------------
    skip_hash = has_double_hash(sumstat_file)
    #print(f"{'Has' if skip_hash else 'No'} '##' metadata")
    # ---- line count ---------------------------------------------------
    shell_cnt = count_data_lines(sumstat_file, skip_hash)
    #print(f"Data lines (count): {shell_cnt:,}")
    # ---- Polars read --------------------------------------------------
    if Path(sumstat_file).suffix.lower() == ".zip":
        try:
            with zipfile.ZipFile(sumstat_file) as zf:
                csvs = [f for f in zf.namelist() if f.lower().endswith(('.csv', '.tsv'))]
                if not csvs:
                    raise ValueError("No CSV/TSV in .zip")
                if len(csvs) > 1:
                    print(f"Multiple files; using: {csvs[0]}")
                txt = zf.read(csvs[0]).decode(errors="ignore")
                df = pl.read_csv(
                    BytesIO(txt.encode("utf-8")),
                    separator=delim,
                    comment_prefix="##" if skip_hash else None,
                    has_header=True,
                    ignore_errors=True,
                    null_values=["NA", "na", ".", ""],
                    truncate_ragged_lines=True,
                )
        except Exception as e:
            raise RuntimeError(f"Failed to read .zip: {e}")
    else:
        try:
            df = pl.read_csv(
                sumstat_file,
                separator=delim,
                comment_prefix="##" if skip_hash else None,
                has_header=True,
                ignore_errors=True,
                null_values=["NA", "na", ".", ""],
                truncate_ragged_lines=True,
            )
        except Exception as e:
            raise RuntimeError(f"Polars failed: {e}")
    #print("Polars read successful")
    # ---- clean --------------------------------------------------------
    df = df.rename({c: c.strip() for c in df.columns})
    for col in df.columns:
        if df[col].dtype in (pl.Utf8, pl.String):
            df = df.with_columns(pl.col(col).str.strip_chars())
    # ---- QC -----------------------------------------------------------
    polars_rows = df.height
    #print(f"Final DataFrame: {polars_rows:,} rows × {len(df.columns)} columns")
    #print(f"QC → Shell count: {shell_cnt:,} | Polars rows: {polars_rows:,}")
    if abs(shell_cnt - polars_rows) > 1:
        print(
            f"⚠️ Warning: Mismatch between shell line-count ({shell_cnt}) and Polars row-count ({polars_rows}).\n"
            "   • This is usually expected when one method counts the header line and the other does not.\n"
            "   • If the difference is large, please verify the input file for formatting issues.\n"
        )
    return df, shell_cnt, polars_rows



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
    print(f"   Tried constructed path: {constructed_path}")
    return "NA"





def read_config(csv_path):
    """
    Convert a CSV configuration file into a JSON configuration file,
    trim whitespace around keys and values, and return as a list of dicts.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV configuration file.
    output_dir : str
        Directory where the JSON configuration file will be saved.

    Returns
    -------
    list[dict]
        List of configuration dictionaries with trimmed keys and values.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"❌ The file '{csv_path}' was not found.")
    output_dir=os.path.dirname(csv_path)
    os.makedirs(output_dir, exist_ok=True)
    # Read the CSV file
    df = pd.read_csv(csv_path)
    df=df.fillna('NA')
    # Strip whitespace from column names and all string cells
    df.columns = df.columns.str.strip()
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Define full JSON output path
    json_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".json"
    json_path = os.path.join(output_dir, json_filename)
    # Convert to JSON and save
    df.to_json(json_path, orient="records", indent=4)
    #print(f"✅ Saved JSON config to: {json_path}")
    # Load JSON and ensure keys and values are stripped (extra safety)
    with open(json_path, "r") as f:
        cfg_list = json.load(f)
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
