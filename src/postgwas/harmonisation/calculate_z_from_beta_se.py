import polars as pl
from pathlib import Path
import io
from typing import Tuple, Dict


import sys

def calculate_z_from_beta_se(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Compute imputed Z-score ('imp_z_col') from BETA and SE columns.

    Formula:
        imp_z_col = BETA / SE

    This version:
      • Silent (no screen printing except summary)
      • Writes full logs to: logs/{gwas_outputname}_chr{chromosome}_calculate_z.log
      • Safe for multiprocessing
    """

    # -------------------------------------------------------
    # Setup logfile
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

    # -------------------------------------------------------
    # Internal log buffer (multiprocessing-safe)
    # -------------------------------------------------------
    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print("\n🧩 Starting computation of imputed Z-scores (imp_z_col)...")

    beta_col = sample_column_dict.get("beta_col", "NA")
    se_col   = sample_column_dict.get("se_col", "NA")

    # -------------------------------------------------------
    # Validate input columns
    # -------------------------------------------------------
    if (
        beta_col == "NA" or
        se_col == "NA" or
        beta_col not in df.columns or
        se_col not in df.columns
    ):
        raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

    # -------------------------------------------------------
    # Cast to numeric
    # -------------------------------------------------------
    df = df.with_columns([
        pl.col(beta_col).cast(pl.Float64, strict=False),
        pl.col(se_col).cast(pl.Float64, strict=False)
    ])

    # -------------------------------------------------------
    # Filtering conditions (with finite checks)
    # -------------------------------------------------------
    cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
    cond_se_null   = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
    cond_beta_zero = (pl.col(beta_col) == 0)
    cond_se_invalid = (pl.col(se_col) <= 0)

    # -------------------------------------------------------
    # Partition removed variants
    # -------------------------------------------------------
    df_beta_null = df.filter(cond_beta_null)
    df_se_null   = df.filter(~cond_beta_null & cond_se_null)
    df_beta_zero = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
    df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

    # -------------------------------------------------------
    # Combine removed variants (safe)
    # -------------------------------------------------------
    df_removed = pl.concat([
        df_beta_null,
        df_se_null,
        df_beta_zero,
        df_se_invalid
    ]).unique()

    # -------------------------------------------------------
    # Keep valid variants
    # -------------------------------------------------------
    condition_valid = (
        (~cond_beta_null) &
        (~cond_se_null) &
        (~cond_beta_zero) &
        (~cond_se_invalid)
    )

    df = df.filter(condition_valid)

    # -------------------------------------------------------
    # Counts
    # -------------------------------------------------------
    n_initial = qc_info["initial_variants"]

    n_beta_null = df_beta_null.height
    n_se_null   = df_se_null.height
    n_beta_zero = df_beta_zero.height
    n_se_invalid = df_se_invalid.height

    n_removed = n_beta_null + n_se_null + n_beta_zero + n_se_invalid
    n_final = df.height

    pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

    # -------------------------------------------------------
    # Print summary
    # -------------------------------------------------------
    indent = "\t" * 4  # 4-level indentation

    msg = (
        f"\n{indent}📊 [Chr {chromosome}] Removed {n_removed:,} variants with invalid Beta/SE.\n"
        f"{indent}{'-'*50}\n"
        f"{indent}• Initial variants : {n_initial:,}\n"
        f"{indent}• Removed total    : {n_removed:,} ({pct_removed:.2f}%)\n"
        f"{indent}    ├─ Null BETA   : {n_beta_null:,}\n"
        f"{indent}    ├─ Null SE     : {n_se_null:,}\n"
        f"{indent}    ├─ BETA = 0    : {n_beta_zero:,}\n"
        f"{indent}    └─ SE <= 0     : {n_se_invalid:,}\n"
        f"{indent}• Remaining        : {n_final:,}\n"
        f"{indent}{'-'*50}\n"
    )
    sys.stdout.write(msg)
    sys.stdout.flush()

    # -------------------------------------------------------
    # Log summary
    # -------------------------------------------------------
    log_print(f"[Chr {chromosome}] QC Summary:")
    log_print(f"Initial: {n_initial}")
    log_print(f"Removed: {n_removed} ({pct_removed:.2f}%)")
    log_print(f"Null BETA: {n_beta_null}")
    log_print(f"Null SE: {n_se_null}")
    log_print(f"BETA=0: {n_beta_zero}")
    log_print(f"SE<=0: {n_se_invalid}")
    log_print(f"Remaining: {n_final}")

    # -------------------------------------------------------
    # QC info update
    # -------------------------------------------------------
    qc_info.update({
        "variants_removed_due_to_null_beta": n_beta_null,
        "variants_removed_due_to_null_se": n_se_null,
        "variants_removed_due_to_zero_beta": n_beta_zero,
        "variants_removed_due_to_invalid_se": n_se_invalid,
        "variants_removed_invalid_beta_se": n_removed,
        "variants_after_filter_invalid_beta_se": n_final
    })

    # -------------------------------------------------------
    # Compute Z safely
    # -------------------------------------------------------
    log_print("📈 Computing imp_z_col = BETA / SE …")

    df = df.with_columns(
        pl.when(pl.col(se_col) > 0)
        .then(pl.col(beta_col) / pl.col(se_col))
        .otherwise(None)
        .alias("imp_z_col")
    )

    sample_column_dict["imp_z_col"] = "imp_z_col"

    # -------------------------------------------------------
    # Z summary
    # -------------------------------------------------------
    z_summary = df.select([
        pl.col("imp_z_col").min().alias("z_min"),
        pl.col("imp_z_col").max().alias("z_max"),
        pl.col("imp_z_col").mean().alias("z_mean"),
        pl.col("imp_z_col").std().alias("z_std"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    log_print(
        f"✅ Imputed Z-score summary: "
        f"min={z_summary['z_min']:.6f}, "
        f"max={z_summary['z_max']:.6f}, "
        f"mean={z_summary['z_mean']:.6f}, "
        f"std={z_summary['z_std']:.6f} "
        f"(n={z_summary['total']:,})"
    )

    qc_info.update({
        "z_min": z_summary["z_min"],
        "z_max": z_summary["z_max"],
        "z_mean": z_summary["z_mean"],
        "z_std": z_summary["z_std"],
        "total_variants": z_summary["total"],
    })

    log_print("🎯 Computation of imp_z_col completed.\n")

    # -------------------------------------------------------
    # Write logfile
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, qc_info, sample_column_dict

# def calculate_z_from_beta_se(
#     chromosome: str,
#     df: pl.DataFrame,
#     sample_column_dict: dict
# ) -> Tuple[pl.DataFrame, Dict, dict]:
#     """
#     Compute imputed Z-score ('imp_z_col') from BETA and SE columns.

#     Formula:
#         imp_z_col = BETA / SE

#     This version:
#       • Silent (no screen printing except summary)
#       • Writes full logs to: logs/{gwas_outputname}_chr{chromosome}_calculate_z.log
#       • Safe for multiprocessing
#     """

#     # -------------------------------------------------------
#     # Setup logfile
#     # -------------------------------------------------------
#     gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
#     output_dir      = sample_column_dict.get("output_folder", ".")
#     log_dir         = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)

#     log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_zscore_from_beta_and_se.log"

#     # -------------------------------------------------------
#     # Internal log buffer (multiprocessing-safe)
#     # -------------------------------------------------------
#     log_buffer = io.StringIO()

#     def log_print(*args):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n")

#     # -------------------------------------------------------
#     # MAIN LOGIC
#     # -------------------------------------------------------
#     qc_info = {"initial_variants": df.height}

#     log_print("\n🧩 Starting computation of imputed Z-scores (imp_z_col)...")

#     beta_col = sample_column_dict.get("beta_col", "NA")
#     se_col   = sample_column_dict.get("se_col", "NA")

#     # -------------------------------------------------------
#     # Validate input columns
#     # -------------------------------------------------------
#     if (
#         beta_col == "NA" or
#         se_col == "NA" or
#         beta_col not in df.columns or
#         se_col not in df.columns
#     ):
#         raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

#     # -------------------------------------------------------
#     # Cast to numeric
#     # -------------------------------------------------------
#     df = df.with_columns([
#         pl.col(beta_col).cast(pl.Float64, strict=False),
#         pl.col(se_col).cast(pl.Float64, strict=False)
#     ])

#     # -------------------------------------------------------
#     # Filtering conditions (with finite checks)
#     # -------------------------------------------------------
#     cond_beta_null = pl.col(beta_col).is_null() | ~pl.col(beta_col).is_finite()
#     cond_se_null   = pl.col(se_col).is_null()   | ~pl.col(se_col).is_finite()
#     cond_beta_zero = (pl.col(beta_col) == 0)
#     cond_se_invalid = (pl.col(se_col) <= 0)

#     # -------------------------------------------------------
#     # Partition removed variants
#     # -------------------------------------------------------
#     df_beta_null = df.filter(cond_beta_null)
#     df_se_null   = df.filter(~cond_beta_null & cond_se_null)
#     df_beta_zero = df.filter(~cond_beta_null & ~cond_se_null & cond_beta_zero)
#     df_se_invalid = df.filter(~cond_beta_null & ~cond_se_null & ~cond_beta_zero & cond_se_invalid)

#     # -------------------------------------------------------
#     # Combine removed variants (safe)
#     # -------------------------------------------------------
#     df_removed = pl.concat([
#         df_beta_null,
#         df_se_null,
#         df_beta_zero,
#         df_se_invalid
#     ]).unique()

#     # -------------------------------------------------------
#     # Keep valid variants
#     # -------------------------------------------------------
#     condition_valid = (
#         (~cond_beta_null) &
#         (~cond_se_null) &
#         (~cond_beta_zero) &
#         (~cond_se_invalid)
#     )

#     df = df.filter(condition_valid)

#     # -------------------------------------------------------
#     # Counts
#     # -------------------------------------------------------
#     n_initial = qc_info["initial_variants"]

#     n_beta_null = df_beta_null.height
#     n_se_null   = df_se_null.height
#     n_beta_zero = df_beta_zero.height
#     n_se_invalid = df_se_invalid.height

#     n_removed = n_beta_null + n_se_null + n_beta_zero + n_se_invalid
#     n_final = df.height

#     pct_removed = (n_removed / n_initial * 100) if n_initial > 0 else 0.0

#     # -------------------------------------------------------
#     # Print summary
#     # -------------------------------------------------------
#     indent = "\t" * 4  # 4-level indentation

#     print(f"\n{indent}📊 [Chr {chromosome}] Removed {n_removed:,} variants with invalid Beta/SE.")
#     print(f"{indent}" + "-"*50)

#     print(f"{indent}• Initial variants : {n_initial:,}")
#     print(f"{indent}• Removed total    : {n_removed:,} ({pct_removed:.2f}%)")

#     print(f"{indent}    ├─ Null BETA   : {n_beta_null:,}")
#     print(f"{indent}    ├─ Null SE     : {n_se_null:,}")
#     print(f"{indent}    ├─ BETA = 0    : {n_beta_zero:,}")
#     print(f"{indent}    └─ SE <= 0     : {n_se_invalid:,}")

#     print(f"{indent}• Remaining        : {n_final:,}")
#     print(f"{indent}" + "-"*50)

#     # -------------------------------------------------------
#     # Log summary
#     # -------------------------------------------------------
#     log_print(f"[Chr {chromosome}] QC Summary:")
#     log_print(f"Initial: {n_initial}")
#     log_print(f"Removed: {n_removed} ({pct_removed:.2f}%)")
#     log_print(f"Null BETA: {n_beta_null}")
#     log_print(f"Null SE: {n_se_null}")
#     log_print(f"BETA=0: {n_beta_zero}")
#     log_print(f"SE<=0: {n_se_invalid}")
#     log_print(f"Remaining: {n_final}")

#     # -------------------------------------------------------
#     # QC info update
#     # -------------------------------------------------------
#     qc_info.update({
#         "variants_removed_due_to_null_beta": n_beta_null,
#         "variants_removed_due_to_null_se": n_se_null,
#         "variants_removed_due_to_zero_beta": n_beta_zero,
#         "variants_removed_due_to_invalid_se": n_se_invalid,
#         "variants_removed_invalid_beta_se": n_removed,
#         "variants_after_filter_invalid_beta_se": n_final
#     })

#     # -------------------------------------------------------
#     # Compute Z safely
#     # -------------------------------------------------------
#     log_print("📈 Computing imp_z_col = BETA / SE …")

#     df = df.with_columns(
#         pl.when(pl.col(se_col) > 0)
#         .then(pl.col(beta_col) / pl.col(se_col))
#         .otherwise(None)
#         .alias("imp_z_col")
#     )

#     sample_column_dict["imp_z_col"] = "imp_z_col"

#     # -------------------------------------------------------
#     # Z summary
#     # -------------------------------------------------------
#     z_summary = df.select([
#         pl.col("imp_z_col").min().alias("z_min"),
#         pl.col("imp_z_col").max().alias("z_max"),
#         pl.col("imp_z_col").mean().alias("z_mean"),
#         pl.col("imp_z_col").std().alias("z_std"),
#         pl.len().alias("total")
#     ]).to_dicts()[0]

#     log_print(
#         f"✅ Imputed Z-score summary: "
#         f"min={z_summary['z_min']:.6f}, "
#         f"max={z_summary['z_max']:.6f}, "
#         f"mean={z_summary['z_mean']:.6f}, "
#         f"std={z_summary['z_std']:.6f} "
#         f"(n={z_summary['total']:,})"
#     )

#     qc_info.update({
#         "z_min": z_summary["z_min"],
#         "z_max": z_summary["z_max"],
#         "z_mean": z_summary["z_mean"],
#         "z_std": z_summary["z_std"],
#         "total_variants": z_summary["total"],
#     })

#     log_print("🎯 Computation of imp_z_col completed.\n")

#     # -------------------------------------------------------
#     # Write logfile
#     # -------------------------------------------------------
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return df, qc_info, sample_column_dict
