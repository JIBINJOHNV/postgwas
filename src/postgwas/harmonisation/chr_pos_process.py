from pathlib import Path
import io
import os
import polars as pl
from typing import Dict, Tuple


def fix_chr_pos_column(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    drop_mt: bool = True
) -> pl.DataFrame:
    """
    Process chromosome and position columns in a GWAS DataFrame.
    """
    # -------------------------------------------------------
    # Setup logfile
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_fix_chr_pos_column.log"

    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------

    log_print("\n🧬 Starting chromosome/position harmonization...")

    chr_col = sample_column_dict.get("chr_col")
    pos_col = sample_column_dict.get("pos_col")
    chr_pos_col = sample_column_dict.get("chr_pos_col")

    # Normalize NA
    if chr_col in [None, "NA"]:
        chr_col = None
    if pos_col in [None, "NA"]:
        pos_col = None

    # -------------------------------------------------------
    # Case 1 — chr_pos combined column
    # -------------------------------------------------------
    if chr_pos_col and chr_pos_col in df.columns and not chr_col and not pos_col:
        log_print("🔹 Splitting combined chromosome-position column...")

        df = df.with_columns(
            pl.col(chr_pos_col)
            .str.replace_all("[:-]", "_")
            .alias(chr_pos_col)
        )

        df = df.with_columns([
            pl.col(chr_pos_col).str.split("_").list.get(0).alias("Chr"),
            pl.col(chr_pos_col).str.split("_").list.get(1).alias("Pos")
        ])

        sample_column_dict["chr_col"] = "Chr"
        sample_column_dict["pos_col"] = "Pos"
        sample_column_dict["chr_pos_col"] =None
        chr_col, pos_col = "Chr", "Pos"

    # -------------------------------------------------------
    # Cast initial types
    # -------------------------------------------------------
    df = df.with_columns([
        pl.col(chr_col).cast(pl.Utf8).alias(chr_col),
        pl.col(pos_col).cast(pl.Int64, strict=False).alias(pos_col)
    ])

    if not chr_col or chr_col not in df.columns:
        raise ValueError("❌ Chromosome column missing after normalization.")

    # -------------------------------------------------------
    # Normalize chromosome labels
    # -------------------------------------------------------
    log_print("🧬 Normalizing chromosome values...")

    df = df.with_columns(
        pl.col(chr_col)
        .str.replace("^chr", "", literal=False)
        .str.strip_chars()
        .str.to_uppercase()
        .alias(chr_col)
    )

    # Convert "1.0" → "1" safely
    df = df.with_columns(
        pl.when(
            pl.col(chr_col).cast(pl.Float64, strict=False).is_not_null()
        )
        .then(
            pl.col(chr_col)
            .cast(pl.Float64, strict=False)
            .cast(pl.Int64, strict=False)
            .cast(pl.String)
        )
        .otherwise(pl.col(chr_col))
        .alias(chr_col)
    )

    # Convert 23→X, 24→Y, 25→MT
    df = df.with_columns(
        pl.when(pl.col(chr_col) == "23").then(pl.lit("X"))
        .when(pl.col(chr_col) == "24").then(pl.lit("Y"))
        .when(pl.col(chr_col) == "25").then(pl.lit("MT"))
        .otherwise(pl.col(chr_col))
        .alias(chr_col)
    )

    # -------------------------------------------------------
    # Filter allowed chromosomes
    # -------------------------------------------------------
    allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]
    if not drop_mt:
        allowed.append("MT")

    df = df.filter(pl.col(chr_col).is_in(allowed))

    # Remove nulls
    if pos_col in df.columns:
        df = df.filter(~pl.col(chr_col).is_null() & ~pl.col(pos_col).is_null())
    else:
        df = df.filter(~pl.col(chr_col).is_null())

    # -------------------------------------------------------
    # Fix position column only (safe)
    # -------------------------------------------------------
    if pos_col in df.columns:
        df = df.with_columns(
            pl.col(pos_col)
            .cast(pl.Float64, strict=False)
            .cast(pl.Int64, strict=False)
            .alias(pos_col)
        )

    log_print(f"✅ Chromosome normalization complete — {df.height:,} rows remain.\n")

    # -------------------------------------------------------
    # Write log
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, sample_column_dict


# def fix_chr_pos_column(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict,
#     drop_mt: bool = True
# ) -> pl.DataFrame:
#     """
#     Process chromosome and position columns in a GWAS DataFrame.
#     Automatically splits chr_pos_col if needed and normalizes chromosome labels.
#     Converts 23→X, 24→Y, optionally removes MT/25.
#     Multiprocessing-safe version:
#         • No stdout printing
#         • All prints go to logs/<gwas>_chr{chromosome}_chrpos.log
#     """
#     # -------------------------------------------------------
#     # Setup logfile
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir = sample_column_dict.get("output_folder", ".")
#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)
#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_fix_chr_pos_column.log"

#     log_buffer = io.StringIO()

#     def log_print(*args):
#         msg = " ".join(str(a) for a in args)
#         # no real print:
#         log_buffer.write(msg + "\n")

#     # -------------------------------------------------------
#     # MAIN LOGIC (all prints → log_print)
#     # -------------------------------------------------------

#     log_print("\n🧬 Starting chromosome/position harmonization...")

#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")
#     chr_pos_col = sample_column_dict.get("chr_pos_col")

#     # Normalize NA
#     if chr_col in [None, "NA"]:
#         chr_col = None
#     if pos_col in [None, "NA"]:
#         pos_col = None

#     # -------------------------------------------------------
#     # Case 1 — chr_pos combined column
#     # -------------------------------------------------------
#     if chr_pos_col and chr_pos_col in df.columns and not chr_col and not pos_col:
#         log_print("🔹 Splitting combined chromosome-position column...")

#         df = df.with_columns(
#             pl.col(chr_pos_col)
#             .str.replace_all("[:-]", "_")
#             .alias(chr_pos_col)
#         )

#         df = df.with_columns([
#             pl.col(chr_pos_col).str.split("_").list.get(0).alias("Chr"),
#             pl.col(chr_pos_col).str.split("_").list.get(1).alias("Pos")
#         ])

#         sample_column_dict["chr_col"] = "Chr"
#         sample_column_dict["pos_col"] = "Pos"
#         chr_col, pos_col = "Chr", "Pos"

#     # -------------------------------------------------------
#     # Cast data types
#     # -------------------------------------------------------
#     df = df.with_columns([
#         pl.col(chr_col).cast(pl.Utf8).alias(chr_col),
#         pl.col(pos_col).cast(pl.Int64).alias(pos_col)
#     ])

#     if not chr_col or chr_col not in df.columns:
#         raise ValueError("❌ Chromosome column missing after normalization.")

#     # -------------------------------------------------------
#     # Normalize chromosome labels
#     # -------------------------------------------------------
#     log_print("🧬 Normalizing chromosome values...")

#     # Strip “chr” prefix
#     df = df.with_columns(pl.col(chr_col).str.replace("chr", "", literal=True))

#     # Convert "1.0" → "1"
#     float_to_str = {f"{i}.0": str(i) for i in range(1, 23)}
#     df = df.with_columns(
#         pl.when(pl.col(chr_col).is_in(list(float_to_str.keys())))
#         .then(pl.col(chr_col).replace(float_to_str))
#         .otherwise(pl.col(chr_col))
#         .alias(chr_col)
#     )

#     # Convert 23→X, 24→Y, 25→MT
#     try:
#         df = df.with_columns(
#             pl.when(pl.col(chr_col) == "23").then(pl.lit("X"))
#             .when(pl.col(chr_col) == "24").then(pl.lit("Y"))
#             .when(pl.col(chr_col) == "25").then(pl.lit("MT"))
#             .otherwise(pl.col(chr_col))
#             .alias(chr_col)
#         )
#     except Exception as e:
#         log_print(f"⚠️ Skipped numeric conversion for X/Y/MT: {e}")

#     # -------------------------------------------------------
#     # Filter allowed chromosomes
#     # -------------------------------------------------------
#     allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]
#     if not drop_mt:
#         allowed.append("MT")

#     df = df.filter(pl.col(chr_col).is_in(allowed))

#     # Remove null pos
#     if pos_col in df.columns:
#         df = df.filter(~pl.col(chr_col).is_null() & ~pl.col(pos_col).is_null())
#     else:
#         df = df.filter(~pl.col(chr_col).is_null())

#     chr_col = sample_column_dict.get("chr_col")
#     pos_col = sample_column_dict.get("pos_col")

#     # Filter to ensure we only touch columns that actually exist in this dataframe
#     cols_to_fix = []

#     # 1. Fix Chromosome: Float(12.0) -> Int(12) -> Str("12")
#     if chr_col in df.columns:
#         cols_to_fix.append(
#             pl.col(chr_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).cast(pl.String)
#         )

#     # 2. Fix Position: Float(1.5e7) -> Int(15000000)
#     if pos_col in df.columns:
#         cols_to_fix.append(
#             pl.col(pos_col).cast(pl.Float64, strict=False).cast(pl.Int64, strict=False)
#         )

#     # 3. Apply all fixes at once
#     if cols_to_fix:
#         df = df.with_columns(cols_to_fix)

#     log_print(f"✅ Chromosome normalization complete — {df.height:,} rows remain.\n")

#     # -------------------------------------------------------
#     # Write log
#     # -------------------------------------------------------
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())
        
#     return df,sample_column_dict


