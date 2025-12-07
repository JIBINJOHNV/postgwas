import polars as pl
import numpy as np
from typing import Tuple, Dict
from pathlib import Path
import io


def calculate_beta_and_se_from_z(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Compute BETA and SE from Z-score (imp_z_col) using zmaf & Neff formulas.
    SAFE FOR MULTIPROCESSING:
        • No console printing
        • All messages logged to per-chromosome logfile
    """

    # -------------------------------------------------------
    # Create logfile
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")

    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_beta_se_from_z.log"

    # -------------------------------------------------------
    # Logging buffer (NO stdout output)
    # -------------------------------------------------------
    log_buffer = io.StringIO()

    def log_print(*args, **kwargs):
        """Write only to file buffer, never to stdout."""
        log_buffer.write(" ".join(str(a) for a in args) + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print("\n🧩 Starting computation of BETA and SE from Z-score...")

    imp_z_col = sample_column_dict.get("imp_z_col", "NA")
    beta_col  = sample_column_dict.get("beta_or_col", "NA")
    se_col    = sample_column_dict.get("se_col", "NA")

    # -------------------------------------------------------
    # 1️⃣ Already present → skip
    # -------------------------------------------------------
    if beta_col in df.columns and se_col in df.columns:
        log_print("✅ BETA+SE columns already present; skipping computation.")

        summary = df.select([
            pl.col(beta_col).min().alias("beta_min"),
            pl.col(beta_col).max().alias("beta_max"),
            pl.col(beta_col).mean().alias("beta_mean"),
            pl.col(beta_col).std().alias("beta_std"),

            pl.col(se_col).min().alias("se_min"),
            pl.col(se_col).max().alias("se_max"),
            pl.col(se_col).mean().alias("se_mean"),
            pl.col(se_col).std().alias("se_std"),

            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update(summary)
        qc_info["status"] = "existing_beta_se"

        # write log
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        return df, qc_info, sample_column_dict

    # -------------------------------------------------------
    # 2️⃣ Check required inputs
    # -------------------------------------------------------
    if imp_z_col == "NA" or imp_z_col not in df.columns:
        log_print("⚠️ No imputed Z-score column found — skipping.")
        qc_info["status"] = "skipped_no_zscore"

        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        return df, qc_info, sample_column_dict

    if "zmaf" not in df.columns or "Neff" not in df.columns:
        raise ValueError("❌ Missing required columns: 'zmaf' and 'Neff' must be present.")

    # -------------------------------------------------------
    # 3️⃣ Ensure numeric types
    # -------------------------------------------------------
    df = df.with_columns([
        pl.col(imp_z_col).cast(pl.Float64, strict=False),
        pl.col("zmaf").cast(pl.Float64, strict=False),
        pl.col("Neff").cast(pl.Float64, strict=False),
    ])

    # -------------------------------------------------------
    # 4️⃣ Z-score summary
    # -------------------------------------------------------
    z_summary = df.select([
        pl.col(imp_z_col).min().alias("z_min"),
        pl.col(imp_z_col).max().alias("z_max"),
        pl.col(imp_z_col).mean().alias("z_mean"),
        pl.col(imp_z_col).std().alias("z_std"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    log_print(
        f"📊 Z summary: min={z_summary['z_min']:.6f}, "
        f"max={z_summary['z_max']:.6f}, "
        f"mean={z_summary['z_mean']:.6f}, std={z_summary['z_std']:.6f}"
    )

    qc_info.update(z_summary)

    # -------------------------------------------------------
    # 5️⃣ Compute BETA
    # -------------------------------------------------------
    if beta_col == "NA" or beta_col not in df.columns:
        log_print("🧮 Computing BETA from Z-score...")

        df = df.with_columns(
            (
                pl.col(imp_z_col) /
                (2 * pl.col("zmaf") * (1 - pl.col("zmaf")) *
                 (pl.col("Neff") + pl.col(imp_z_col)**2)).sqrt()
            ).alias("BETA")
        )

        beta_col = "BETA"
        sample_column_dict["beta_or_col"] = "BETA"
        qc_info["beta_computed"] = True

    else:
        log_print("ℹ️ BETA already present — skipping.")
        qc_info["beta_computed"] = False

    # -------------------------------------------------------
    # 6️⃣ Compute SE
    # -------------------------------------------------------
    if se_col == "NA" or se_col not in df.columns:
        log_print("🧮 Computing SE from Z-score...")

        df = df.with_columns(
            (
                1 /
                (2 * pl.col("zmaf") * (1 - pl.col("zmaf")) *
                 (pl.col("Neff") + pl.col(imp_z_col)**2)).sqrt()
            ).alias("SE")
        )

        se_col = "SE"
        sample_column_dict["se_col"] = "SE"
        qc_info["se_computed"] = True

    else:
        log_print("ℹ️ SE already present — skipping.")
        qc_info["se_computed"] = False

    # -------------------------------------------------------
    # 7️⃣ Final BETA/SE summary
    # -------------------------------------------------------
    summary = df.select([
        pl.col(beta_col).min().alias("beta_min"),
        pl.col(beta_col).max().alias("beta_max"),
        pl.col(beta_col).mean().alias("beta_mean"),
        pl.col(beta_col).std().alias("beta_std"),

        pl.col(se_col).min().alias("se_min"),
        pl.col(se_col).max().alias("se_max"),
        pl.col(se_col).mean().alias("se_mean"),
        pl.col(se_col).std().alias("se_std"),

        pl.len().alias("total")
    ]).to_dicts()[0]

    qc_info.update(summary)

    log_print(
        f"📈 BETA summary: min={summary['beta_min']:.6e}, "
        f"max={summary['beta_max']:.6e}, mean={summary['beta_mean']:.6e}, std={summary['beta_std']:.6e}"
    )
    log_print(
        f"📉 SE summary:   min={summary['se_min']:.6e}, "
        f"max={summary['se_max']:.6e}, mean={summary['se_mean']:.6e}, std={summary['se_std']:.6e}"
    )

    log_print("🎯 Completed BETA/SE computation from Z-score.\n")

    # -------------------------------------------------------
    # Save logs
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, qc_info, sample_column_dict
