import polars as pl
import io, os
from pathlib import Path
from typing import Tuple, Dict


def add_or_calculate_info(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    info_file: str = "NA",
    info_column: str = "info",
    default_info_file: str = "NA",
    default_info_column: str = "info"
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Add or merge imputation INFO score column into GWAS summary dataframe using Polars.
    Silent version:
        • No printing to screen
        • All logs written to logs/{gwas_outputname}_chr{chromosome}_info.log
        • Fully multiprocessing-safe
    """

    # -------------------------------------------------------
    # Setup chromosome-specific log file
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir = sample_column_dict.get("output_folder", ".")
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_add_or_calculate_info.log"

    # Multiprocessing-safe, SILENT buffered logger
    log_buffer = io.StringIO()

    def log_print(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")   # <-- NO printing to screen

    # -------------------------------------------------------
    # MAIN LOGIC
    # -------------------------------------------------------
    qc_info = {"initial_variants": df.height}

    log_print("\n🧩 Starting INFO score harmonization...")

    # Extract mapping
    chr_col = sample_column_dict["chr_col"]
    pos_col = sample_column_dict["pos_col"]
    ea_col = sample_column_dict["ea_col"]
    oa_col = sample_column_dict["oa_col"]
    imp_info_col = sample_column_dict.get("imp_info_col", "NA")

    # ===============================================================
    # STEP 1: Use internal INFO column if present
    # ===============================================================
    if imp_info_col != "NA" and imp_info_col in df.columns:
        log_print(f"🔍 Found internal INFO column '{imp_info_col}' — validating numeric values.")

        df = df.with_columns(pl.col(imp_info_col).cast(pl.Float64, strict=False))

        invalid_before = df.filter(
            (pl.col(imp_info_col) < 0.0) | (pl.col(imp_info_col) > 1.05)
        ).height

        df = df.with_columns(pl.col(imp_info_col).clip(0.0, 1.0).alias(imp_info_col))

        stats = df.select([
            pl.col(imp_info_col).min().alias("min"),
            pl.col(imp_info_col).max().alias("max"),
            pl.col(imp_info_col).mean().alias("mean"),
            pl.col(imp_info_col).median().alias("median"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        qc_info.update({
            "info_source": "internal",
            "invalid_values_clipped": invalid_before,
            "min_info": stats["min"],
            "max_info": stats["max"],
            "mean_info": stats["mean"],
            "median_info": stats["median"],
            "variants_before": qc_info["initial_variants"],
            "variants_after": stats["total"],
        })

        log_print(
            f"✅ INFO (internal): range={stats['min']:.3f}–{stats['max']:.3f}, "
            f"mean={stats['mean']:.3f}, median={stats['median']:.3f}, n={stats['total']:,}"
        )

        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())

        return df, qc_info, sample_column_dict

    # ===============================================================
    # STEP 2: Choose external INFO source
    # ===============================================================
    if info_file != "NA" and os.path.exists(info_file):
        info_file_to_use = info_file
        info_col_to_use = info_column
        info_source = f"user provided"
    elif default_info_file != "NA" and os.path.exists(default_info_file):
        info_file_to_use = default_info_file
        info_col_to_use = default_info_column
        info_source = f"default"
    else:
        raise FileNotFoundError(
            "❌ No valid INFO source provided (neither internal INFO nor external INFO files)."
        )

    log_print(
        f"📂 Using {info_source} INFO file: {os.path.basename(info_file_to_use)} "
        f"(column: '{info_col_to_use}')"
    )

    # ===============================================================
    # Helper: Load external INFO file
    # ===============================================================
    def load_info_file(path: str, colname: str) -> pl.DataFrame:
        info_df = (
            pl.read_csv(path, separator=r"\s+")
            .select(["chr", "pos", "a1", "a2", colname])
            .with_columns([
                pl.col("pos").cast(pl.Int64),
                pl.col("chr").cast(pl.Utf8).str.replace("0X", "X"),
                pl.col("a1").str.to_uppercase(),
                pl.col("a2").str.to_uppercase(),
                pl.col(colname).cast(pl.Float64, strict=False)
            ])
            .rename({"chr": chr_col, "pos": pos_col, "a1": ea_col, "a2": oa_col})
            .with_columns(pl.col(colname).clip(0.0, 1.0).alias(colname))
        )
        log_print(f"✅ INFO file loaded → {info_df.height:,} variants.")
        return info_df

    # ===============================================================
    # STEP 3: Merge INFO (direct + flipped)
    # ===============================================================
    info_df = load_info_file(info_file_to_use, info_col_to_use)

    # Normalize alleles in df
    df = df.with_columns([
        pl.col(ea_col).cast(pl.Utf8).str.to_uppercase(),
        pl.col(oa_col).cast(pl.Utf8).str.to_uppercase(),
    ])

    # Direct match
    df1 = df.join(info_df, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

    # Flipped match
    info_df_flip = (
        info_df.rename({ea_col: oa_col, oa_col: ea_col})
        .with_columns((1.0 - pl.col(info_col_to_use)).alias(info_col_to_use))
    )

    df2 = df.join(info_df_flip, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

    # Combine & deduplicate
    df_merged = pl.concat([df1, df2], how="vertical").unique(
        subset=[chr_col, pos_col, ea_col, oa_col]
    )

    # ===============================================================
    # STEP 4: Statistics
    # ===============================================================
    stats = df_merged.select([
        pl.col(info_col_to_use).min().alias("min"),
        pl.col(info_col_to_use).max().alias("max"),
        pl.col(info_col_to_use).mean().alias("mean"),
        pl.col(info_col_to_use).median().alias("median"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    qc_info.update({
        "info_source": info_source,
        "info_file_used": info_file_to_use,
        "info_column_used": info_col_to_use,
        "min_info": stats["min"],
        "max_info": stats["max"],
        "mean_info": stats["mean"],
        "median_info": stats["median"],
        "variants_before": qc_info["initial_variants"],
        "variants_after": stats["total"],
    })

    sample_column_dict["imp_info_col"] = info_col_to_use

    log_print(
        f"✅ INFO merged ({info_source}): {stats['total']:,} variants "
        f"({stats['min']:.3f}–{stats['max']:.3f}, mean={stats['mean']:.3f}, median={stats['median']:.3f})"
    )
    log_print(
        f"📊 Variants before={qc_info['initial_variants']:,}, after={stats['total']:,}"
    )
    log_print("🎯 INFO harmonization complete.\n")

    # -------------------------------------------------------
    # Write log buffer to log file
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df_merged, qc_info, sample_column_dict
