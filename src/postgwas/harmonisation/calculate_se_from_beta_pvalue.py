import polars as pl
from scipy.stats import norm
from typing import Tuple, Dict
import numpy as np
import polars as pl
from typing import Tuple, Dict
from pathlib import Path
import numpy as np
from scipy.stats import norm
import io





def calculate_se_from_beta_pvalue(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    tail: int = 2
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Compute Standard Error (SE) from BETA and processed P-values when SE is missing.
    This version includes dedicated multiprocessing-safe per-chromosome logging.
    """

    # -------------------------------------------------------
    # Create logfile path
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_calculate_se_from_beta_and_pvalue.log"

    # -------------------------------------------------------
    # INTERNAL LOG BUFFER (safe in multiprocessing)
    # -------------------------------------------------------
    log_buffer = io.StringIO()

    def log_print(*args):
        """Write only to buffer — never to console."""
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print("\n🧩 Starting SE calculation from BETA and processed p-values...")

    beta_col = sample_column_dict.get("beta_or_col", "NA")
    se_col   = sample_column_dict.get("se_col", "NA")
    pval_col = sample_column_dict.get("pval_col", "NA")

    # --- Validate required columns ---
    if (
        beta_col == "NA" or
        pval_col == "NA" or
        beta_col not in df.columns or
        pval_col not in df.columns
    ):
        raise ValueError("❌ Missing required columns: BETA and/or p-value for SE calculation.")

    # === Case: SE is already provided ===
    if se_col != "NA" and se_col in df.columns:
        log_print(f"ℹ️ SE column '{se_col}' already exists — skipping recalculation.")
        qc_info["status"] = "SE already present"

        stats = df.select([
            pl.col(se_col).min().alias("min"),
            pl.col(se_col).max().alias("max"),
            pl.col(se_col).mean().alias("mean"),
            pl.col(se_col).median().alias("median"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update({
            "se_min": stats["min"],
            "se_max": stats["max"],
            "se_mean": stats["mean"],
            "se_median": stats["median"],
            "total_variants": stats["total"],
        })

        # WRITE LOG FILE
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        return df, qc_info, sample_column_dict

    # --- Cast required columns safely ---
    df = df.with_columns([
        pl.col(beta_col).cast(pl.Float64, strict=False),
        pl.col(pval_col).cast(pl.Float64, strict=False)
    ])

    log_print("📈 Computing SE using scipy.stats.norm.ppf()...")

    # Convert to pandas for norm.ppf
    pdf = df.select([beta_col, pval_col]).to_pandas()

    # Protect norm.ppf against infinite results
    pdf[pval_col] = pdf[pval_col].clip(lower=1e-300, upper=0.999999999)

    # Main SE formula
    pdf["SE"] = np.abs(pdf[beta_col] / norm.ppf(pdf[pval_col] / tail))

    # Merge back into Polars
    df = df.with_columns(pl.Series("SE", pdf["SE"].astype(float)))
    sample_column_dict["se_col"] = "SE"

    # --- QC summary ---
    stats = df.select([
        pl.col("SE").min().alias("min"),
        pl.col("SE").max().alias("max"),
        pl.col("SE").mean().alias("mean"),
        pl.col("SE").median().alias("median"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    qc_info.update({
        "se_min": stats["min"],
        "se_max": stats["max"],
        "se_mean": stats["mean"],
        "se_median": stats["median"],
        "total_variants": stats["total"],
        "status": "SE calculated from BETA and p-value"
    })

    log_print(
        f"✅ SE calculated: min={stats['min']:.6f}, max={stats['max']:.6f}, "
        f"mean={stats['mean']:.6f}, median={stats['median']:.6f} "
        f"over {stats['total']:,} variants."
    )

    log_print("🎯 SE calculation complete.\n")

    # -------------------------------------------------------
    # WRITE LOG FILE
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df, qc_info, sample_column_dict
