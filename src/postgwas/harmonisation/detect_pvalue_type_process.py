import polars as pl
import math
from typing import Dict, Tuple, Any,Union
import io
from pathlib import Path


def detect_pval_type(
    df: pl.DataFrame,
    sample_column_dict,
    proportion_threshold: float = 0.10):
    
    pval_col = sample_column_dict.get("pval_col")
    if not pval_col or pval_col not in df.columns:
        raise ValueError("Missing or invalid p-value column in sample_column_dict.")
    pvals = (
        df.select(pl.col(pval_col).cast(pl.Float64, strict=False))
        .drop_nulls()
        .rename({pval_col: "p"})
    )
    total = pvals.height
    if total == 0:
        raise ValueError("No numeric P-values found for detection.")
    over_1p2 = pvals.filter(pl.col("p") > 1.2).height / total
    within_0_1 = pvals.filter((pl.col("p") >= 0) & (pl.col("p") <= 1.05)).height / total
    max_val = float(pvals.select(pl.col("p").max()).item())
    median_val = float(pvals.select(pl.col("p").median()).item())
    
    if over_1p2 >= proportion_threshold:
        detected_type = "mlogp"
    elif within_0_1 >= 0.95:
        detected_type = "pvalue"
    elif max_val > 10:
        detected_type = "mlogp"
    else:
        detected_type = "pvalue"
    return {
        "detected_type": detected_type,
        "n_values": total,
        "over_1.2_fraction": over_1p2,
        "within_[0,1.05]_fraction": within_0_1,
        "median": median_val,
        "max": max_val,
    }



def convert_pval_to_mlogp(
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "LP",
    min_p: float = 1e-300,
    max_p: float = 0.99):
    """
    Safely convert raw p-values (0–1) → -log10(p).
    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing p-values.
    sample_column_dict : dict[str, str]
        Must include 'pval_col' key with the column name.
    output_col : str, optional
        Name for the new -log10(p) column.
    min_p, max_p : float, optional
        Clipping bounds for out-of-range p-values.
    Returns
    -------
    (DataFrame, dict)
        Cleaned DataFrame and summary stats.
    """
    pval_col = sample_column_dict["pval_col"]
    if pval_col not in df.columns:
        raise ValueError(f"Column '{pval_col}' not found in DataFrame.")
    
    # 2️⃣ Remove invalid or missing values
    df = df.filter(
        (pl.col(pval_col).is_not_null()) &
        (pl.col(pval_col) > 0) &
        (pl.col(pval_col) < 1) )
    # 3️⃣ Clip values within bounds
    df = df.with_columns(
        pl.col(pval_col)
          .clip(min_p, max_p)
          .alias("__p_clipped") )
    # 4️⃣ Convert to -log10(p)
    df = df.with_columns(
        (-pl.col("__p_clipped").log10()).clip(lower_bound=0.0).alias(output_col)
    ).drop("__p_clipped")
    
    # 5️⃣ Summary stats
    stats = df.select(
        pl.col(output_col).min().alias("min_LP"),
        pl.col(output_col).max().alias("max_LP"),
        pl.col(output_col).mean().alias("mean_LP"),
        pl.col(output_col).median().alias("median_LP"),
    ).to_dicts()[0]
    print(f"✅ Converted raw p-values → '{output_col}'")
    print(
        f"LP summary: min={stats['min_LP']:.3f}, max={stats['max_LP']:.3f}, "
        f"mean={stats['mean_LP']:.3f}, median={stats['median_LP']:.3f}"
    )
    sample_column_dict["pval_col"] = output_col
    return df, stats,sample_column_dict



def convert_mlogp_to_pval(
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "PVAL",
    min_lp: float = 0.0,
    max_lp: float = 300.0,
):
    """
    Safely convert -log10(p) → raw p-values.
    
    Parameters
    ----------
    df : pl.DataFrame
        Input DataFrame containing -log10(p) values.
    sample_column_dict : dict[str, str]
        Must include 'pval_col' key pointing to the -log10(p) column.
    output_col : str, optional
        Name for the new raw p-value column (default: "PVAL").
    min_lp, max_lp : float, optional
        Bounds for clipping mlogp values before conversion (default: 0–300).
        Prevents overflow for extremely large LP values.
    
    Returns
    -------
    (DataFrame, dict, sample_column_dict)
        Converted DataFrame, summary statistics, and updated sample_column_dict.
    """
    pval_col = sample_column_dict["pval_col"]
    if pval_col not in df.columns:
        raise ValueError(f"Column '{pval_col}' not found in DataFrame.")
    
    # 1️⃣ Filter valid entries
    df = df.filter(pl.col(pval_col).is_not_null())

    # 2️⃣ Clip extreme LP values to avoid underflow
    df = df.with_columns(
        pl.col(pval_col)
          .clip(lower_bound=min_lp, upper_bound=max_lp)
          .alias("__lp_clipped") )
    
    # 3️⃣ Convert -log10(p) → p = 10^(-LP)
    df = df.with_columns(
        (10 ** (-pl.col("__lp_clipped"))).alias(output_col)
    ).drop("__lp_clipped")
    
    # 4️⃣ Compute summary stats
    stats = df.select(
        pl.col(output_col).min().alias("min_PVAL"),
        pl.col(output_col).max().alias("max_PVAL"),
        pl.col(output_col).mean().alias("mean_PVAL"),
        pl.col(output_col).median().alias("median_PVAL"),
    ).to_dicts()[0]
    print(f"✅ Converted -log10(p) → '{output_col}'")
    print(
        f"PVAL summary: min={stats['min_PVAL']:.3e}, max={stats['max_PVAL']:.3e}, "
        f"mean={stats['mean_PVAL']:.3e}, median={stats['median_PVAL']:.3e}"
    )
    # 5️⃣ Update column dictionary
    sample_column_dict["pval_col"] = output_col
    return df, stats, sample_column_dict



import polars as pl
from pathlib import Path
import io
from typing import Tuple, Dict


def detect_and_convert_pval(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict,
    output_col: str = "LP",
    proportion_threshold: float = 0.10,
    max_allowed_lp: float = 300.0,
):
    """
    Logged wrapper around your existing detect_and_convert_pval, which:
        • creates a log file per chromosome
        • captures all prints safely (no stdout hijacking)
        • preserves ALL your commented code
    """

    # -------------------------------------------------------
    # Setup path + log buffer
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    logfile = log_dir / f"{gwas_outputname}_chr{chromosome}_detect_pval.log"

    buffer = io.StringIO()

    def log_print(*args):
        """Write ONLY to log buffer — never to screen."""
        msg = " ".join(str(x) for x in args)
        buffer.write(msg + "\n")

    # -------------------------------------------------------
    # Begin logic (your code preserved)
    # -------------------------------------------------------
    log_print("\n🧩 Starting P-value detection + harmonization…")

    # 1️⃣ Convert to float
    pval_col = sample_column_dict["pval_col"]
    df = df.with_columns(pl.col(pval_col).cast(pl.Float64))

    # detect raw vs -log10(p)
    detection_info = detect_pval_type(df, sample_column_dict, proportion_threshold)
    qc_info = detection_info.copy()

    print(f"        ✅ P value Detected as {detection_info['detected_type']} for chromosome {chromosome}")
    # --- Case: -log10(p) detected ---
    if detection_info["detected_type"] == "mlogp":
        log_print("Converting -log10(p) P-values to raw P-values…")
        df, stats, sample_column_dict = convert_pval_to_mlogp(
            df=df,
            sample_column_dict=sample_column_dict,
            output_col=output_col
        )

        qc_info.update(stats)
        qc_info["conversion"] = "convert_pval_to_mlogp applied"

    # --- Case: already raw p-values ---
    else:
        df = df.with_columns(pl.col(sample_column_dict["pval_col"]).alias(output_col))
        qc_info["conversion"] = "none_needed"
        log_print(f"P-values already raw; stored in '{output_col}'.")

    # ------------------------------------------------------------
    #  PRESERVED — DO NOT REMOVE ANY COMMENTS
    # ------------------------------------------------------------
    # if detection_info["detected_type"] == "pvalue":
    #     print("Converting raw P-values to -log10(p)…")
    #     df, stats,sample_column_dict = convert_pval_to_mlogp(df, sample_column_dict, output_col)
    #     qc_info.update(stats)
    #     qc_info["conversion"] = "raw_to_mlogp_applied"
    # else:
    #     # Already in -log10(p)
    #     df = df.with_columns(pl.col(sample_column_dict["pval_col"]).alias(output_col))
    #     qc_info["conversion"] = "none_needed"
    #     print(f"P-values already in -log10(p); stored in '{output_col}'.")
    # ------------------------------------------------------------

    # ------------------------------------------------------------------
    # Enforce upper bound on -log10(p)
    # ------------------------------------------------------------------
    n_before = df.height

    df = df.with_columns(
        pl.col(output_col).clip(upper_bound=max_allowed_lp).alias(output_col)
    )

    n_clipped = df.select(
        (pl.col(output_col) == max_allowed_lp).sum()
    ).item()

    if n_clipped:
        log_print(
            f"Clipped {n_clipped:,}/{n_before:,} {output_col} values to "
            f"{max_allowed_lp:.1f}  (p ≤ 10⁻{max_allowed_lp:.0f})"
        )

    qc_info["max_lp_clipped"] = max_allowed_lp
    qc_info["n_lp_clipped"]   = n_clipped

    sample_column_dict["pval_col"] = output_col

    log_print("P-value harmonisation complete.\n")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    stats = df.select([
        pl.col(pval_col).min().alias("originalPValue_min"),
        pl.col(pval_col).max().alias("originalPValue_max"),
        pl.col(pval_col).mean().alias("originalPValue_mean"),
        pl.col(pval_col).median().alias("originalPValue_median"),
        pl.len().alias("originalPValue_total"),

        pl.col(output_col).min().alias("Formated_PValue_min"),
        pl.col(output_col).max().alias("Formated_PValue_max"),
        pl.col(output_col).mean().alias("Formated_PValue_mean"),
        pl.col(output_col).median().alias("Formated_PValue_median"),
        pl.len().alias("Formated_PValue_total"),
    ]).to_dicts()[0]

    qc_info.update(stats)
    qc_info["status"] = "P-value formatting summary"

    # -------------------------------------------------------
    # Write logfile
    # -------------------------------------------------------
    with open(logfile, "w") as f:
        f.write(buffer.getvalue())

    return df, qc_info, sample_column_dict
