import polars as pl
from pathlib import Path
import io
from typing import Tuple, Dict

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
      • Silent (no screen printing)
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
        """Write only to log buffer (never to console)."""
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print("\n🧩 Starting computation of imputed Z-scores (imp_z_col)...")

    beta_col = sample_column_dict.get("beta_col", "NA")
    se_col   = sample_column_dict.get("se_col", "NA")


    # --- Validate input columns ---
    if (
        beta_col == "NA" or
        se_col == "NA" or
        beta_col not in df.columns or
        se_col not in df.columns
    ):
        raise ValueError("❌ Missing required columns for Z-score computation: BETA and/or SE.")

    # 1️⃣ Cast first to ensure numeric operations work
    df = df.with_columns([
        pl.col(beta_col).cast(pl.Float64, strict=False),
        pl.col(se_col).cast(pl.Float64, strict=False)
    ])

    # -------------------------------------------------------
    # 1️⃣ Filtering: Remove Nulls, Beta=0, and SE <= 0
    # -------------------------------------------------------
    n_initial = df.height
    
    # Define a filter mask for "invalid" variants
    # We want to keep: Beta != 0 AND SE > 0 AND neither is Null
    df = df.filter(
        (pl.col(beta_col).is_not_null()) &
        (pl.col(se_col).is_not_null()) &
        (pl.col(beta_col) != 0) &
        (pl.col(se_col) > 0)
    )
    
    n_final = df.height
    n_removed = n_initial - n_final

    # Print to screen (as requested)
    print(f"  [Chr {chromosome}] Removed {n_removed:,} variants with invalid Beta/SE (Null, Beta=0, or SE<=0).")
    
    # Log to file
    filter_msg = (
        f"📊 [Chr {chromosome}] Filtering stats for variants with "
        f"Missing Beta/SE, Zero Beta, or Non-positive SE:\n"
        f"   • Initial: {n_initial:,}\n"
        f"   • Removed: {n_removed:,}\n"
        f"   • Remaining: {n_final:,}"
    )

    print(filter_msg)      # Prints to your terminal screen
    log_print(filter_msg) # Writes to your per-chromosome log file

    # Update QC info with your specific keys
    qc_info["variants_removed_invalid_beta_se"] = n_removed
    qc_info["variants_after_filter_invalid_beta_se"] = n_final

    # --- Compute Z ---
    log_print("📈 Computing imp_z_col = BETA / SE …")
    df = df.with_columns((pl.col(beta_col) / pl.col(se_col)).alias("imp_z_col"))
    sample_column_dict["imp_z_col"] = "imp_z_col"

    # --- Summary statistics ---
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
    # WRITE LOGFILE
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, qc_info, sample_column_dict
