import os
import io
import polars as pl
from pathlib import Path
from typing import Tuple, Dict

import os
import io
from pathlib import Path
from typing import Tuple, Dict
import polars as pl


def add_or_calculate_eaf(
    chromosome: str,
    df: pl.DataFrame,
    sample_column_dict: dict,
    eaffile: str = "NA",
    default_eaf_file: str = "NA",
    default_eaf_eafcolumn: str = "EAF",
) -> Tuple[pl.DataFrame, Dict, dict]:
    """
    Multiprocessing-safe version of allele frequency harmonisation.
    All output is written ONLY to:

        logs/{gwas_outputname}_chr{chromosome}_eaf.log

    No output is printed to screen.
    """

    # -------------------------------------------------------
    # Create per-chromosome logfile
    # -------------------------------------------------------
    gwas_outputname = sample_column_dict.get("gwas_outputname", "GWAS")
    output_dir      = sample_column_dict.get("output_folder", ".")
    log_dir         = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_eaf.log"

    # StringIO buffer to collect prints
    log_buffer = io.StringIO()

    def log_print(*args, **kwargs):
        """Write ONLY to logfile, not to stdout."""
        log_buffer.write(" ".join(str(a) for a in args) + "\n")

    # -------------------------------------------------------
    # ORIGINAL LOGIC (unchanged except print → log_print)
    # -------------------------------------------------------

    qc_info = {"initial_variants": df.height}
    log_print("\n🧩 Starting allele frequency harmonization...")

    chr_col = sample_column_dict["chr_col"]
    pos_col = sample_column_dict["pos_col"]
    ea_col  = sample_column_dict["ea_col"]
    oa_col  = sample_column_dict["oa_col"]
    eaf_col = sample_column_dict.get("eaf_col", "NA")
    eafcolumn = sample_column_dict.get("eafcolumn", "EAF")

    # -----------------------------------------------------------------------
    # Helper: Validate EAF column
    # -----------------------------------------------------------------------
    def validate_eaf_column(df_check: pl.DataFrame, colname: str, label: str):
        if colname not in df_check.columns:
            raise ValueError(f"❌ {label} '{colname}' not found in dataframe.")

        total = df_check.height
        out_of_range = df_check.filter(
            (pl.col(colname) < 0.0) | (pl.col(colname) > 1.05)
        ).height
        out_of_range_pct = out_of_range / total if total else 0

        low_freq = df_check.filter(pl.col(colname) <= 0.5).height
        low_freq_pct = low_freq / total if total else 0

        stats = df_check.select([
            pl.col(colname).min().alias("min"),
            pl.col(colname).max().alias("max"),
            pl.col(colname).mean().alias("mean"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        log_print(f"📊 {label} '{colname}' summary before clipping:")
        log_print(f"   Total variants: {stats['total']:,}")
        log_print(f"   Invalid values: {out_of_range:,} ({out_of_range_pct*100:.2f}%)")
        log_print(f"   Range: {stats['min']} – {stats['max']}, Mean: {stats['mean']}")

        if out_of_range_pct > 0.05:
            raise ValueError(
                f"❌ {label} '{colname}' has {out_of_range_pct*100:.2f}% invalid values."
            )
        if low_freq_pct > 0.75:
            raise ValueError(
                f"⚠️ {label} '{colname}' looks like MAF ({low_freq_pct*100:.2f}% ≤ 0.5)."
            )

        df_checked = df_check.with_columns(
            pl.col(colname).clip(0.000001, 1.0).alias(colname)
        )

        clipped_stats = df_checked.select([
            pl.col(colname).min().alias("min"),
            pl.col(colname).max().alias("max"),
            pl.col(colname).mean().alias("mean"),
            pl.len().alias("total")
        ]).to_dicts()[0]

        log_print(f"✅ {label} '{colname}' clipped to [0.000001–1.0]")
        log_print(
            f"   New Range: {clipped_stats['min']} – {clipped_stats['max']}, "
            f"Mean: {clipped_stats['mean']}"
        )

        return df_checked, {
            "out_of_range_pct": out_of_range_pct,
            "low_freq_pct": low_freq_pct,
            "min": clipped_stats["min"],
            "max": clipped_stats["max"],
            "mean": clipped_stats["mean"],
        }

    # -----------------------------------------------------------------------
    # Helper: Load external EAF file
    # -----------------------------------------------------------------------
    def load_eaf_file(path: str) -> pl.DataFrame:
        if path == "NA" or not os.path.exists(path):
            raise FileNotFoundError(f"❌ EAF reference file '{path}' not found.")

        log_print(f"📂 Loading external EAF: {os.path.basename(path)} [expected column: {eafcolumn}]")

        eaf_df = (
            pl.read_csv(path, separator=r"\s+")
            .select(["chr", "pos", "a1", "a2", eafcolumn])
            .with_columns([
                pl.col("pos").cast(pl.Int64),
                pl.col("chr").cast(pl.Utf8).str.replace("0X", "X"),
                pl.col("a1").str.to_uppercase(),
                pl.col("a2").str.to_uppercase(),
            ])
            .rename({"chr": chr_col, "pos": pos_col, "a1": ea_col, "a2": oa_col})
        )

        eaf_df, _ = validate_eaf_column(eaf_df, eafcolumn, "External EAF")
        return eaf_df

    # -----------------------------------------------------------------------
    # STEP 1 — Internal EAF
    # -----------------------------------------------------------------------
    eaf_valid = False

    if eaf_col != "NA" and eaf_col in df.columns:
        log_print(f"🔍 Found internal EAF column '{eaf_col}', validating...")

        try:
            df, result = validate_eaf_column(df, eaf_col, "Internal EAF")
            eaf_valid = True
            qc_info.update({"internal_eaf_valid": True, **result})
        except ValueError as e:
            log_print(f"⚠️ {e}")
            qc_info.update({"internal_eaf_valid": False, "reason": str(e)})
    else:
        log_print("ℹ️ No internal EAF column found.")
        qc_info["internal_eaf_valid"] = False

    # -----------------------------------------------------------------------
    # STEP 2 — If internal valid, compute zMAF
    # -----------------------------------------------------------------------
    if eaf_valid:
        log_print(f"✅ '{eaf_col}' is valid. Computing zMAF directly.")

        df_checked = df.with_columns(
            (pl.when(pl.col(eaf_col) <= 0.5)
             .then(pl.col(eaf_col))
             .otherwise(1 - pl.col(eaf_col)))
            .alias("zmaf")
        )
        sample_column_dict["eaf_col"] = eaf_col

    else:
        # -------------------------------------------------------------------
        # STEP 3 — Load external EAF
        # -------------------------------------------------------------------
        eaf_file_to_use = (
            eaffile if (eaffile != "NA" and os.path.exists(eaffile))
            else default_eaf_file
        )

        if eaf_file_to_use == default_eaf_file:
            eafcolumn = default_eaf_eafcolumn
        else:
            eafcolumn = sample_column_dict.get("eafcolumn", "EAF")

        log_print(f"ℹ️ Using external EAF file: {eaf_file_to_use}")
        log_print(f"📊 Expected EAF column in file: '{eafcolumn}'")

        eaf_df = load_eaf_file(eaf_file_to_use)

        df = df.with_columns([
            pl.col(ea_col).cast(pl.Utf8).str.to_uppercase(),
            pl.col(oa_col).cast(pl.Utf8).str.to_uppercase(),
        ])

        df1 = df.join(eaf_df, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

        eaf_df_flip = (
            eaf_df.rename({ea_col: oa_col, oa_col: ea_col})
            .with_columns((1 - pl.col(eafcolumn)).alias(eafcolumn))
        )

        df2 = df.join(eaf_df_flip, on=[chr_col, pos_col, ea_col, oa_col], how="inner")

        df_checked = pl.concat([df1, df2], how="vertical") \
                       .unique([chr_col, pos_col, ea_col, oa_col])

        df_checked = df_checked.with_columns(
            (pl.when(pl.col(eafcolumn) <= 0.5)
             .then(pl.col(eafcolumn))
             .otherwise(1 - pl.col(eafcolumn)))
            .alias("zmaf")
        )

        sample_column_dict["eaf_col"] = eafcolumn
        qc_info["external_eaf_file"] = eaf_file_to_use

    # -----------------------------------------------------------------------
    # STEP 4 — Final sanity check
    # -----------------------------------------------------------------------
    log_print("🧪 Performing final EAF sanity checks (before and after clipping)...")

    eaf_col_to_check = sample_column_dict["eaf_col"]
    total = df_checked.height

    out_of_bounds_before = df_checked.filter(
        (pl.col(eaf_col_to_check) < 0.000001)
        | (pl.col(eaf_col_to_check) > 1.0)
    ).height

    frac_before = out_of_bounds_before / total if total else 0

    log_print(
        f"⚙️ Before clipping: {out_of_bounds_before:,} "
        f"({frac_before*100:.3f}%) out-of-range."
    )

    df_checked = df_checked.with_columns(
        pl.col(eaf_col_to_check).clip(0.000001, 1.0).alias(eaf_col_to_check)
    )

    stats = df_checked.select([
        pl.col(eaf_col_to_check).min().alias("min"),
        pl.col(eaf_col_to_check).max().alias("max"),
        pl.col(eaf_col_to_check).mean().alias("mean"),
        pl.len().alias("total")
    ]).to_dicts()[0]

    log_print(
        f"✅ After clipping: range={stats['min']} – {stats['max']}, "
        f"mean={stats['mean']} over {stats['total']:,} variants."
    )

    qc_info.update({
        "preclip_out_of_range_count": out_of_bounds_before,
        "preclip_out_of_range_fraction": frac_before,
        "final_min_eaf": stats["min"],
        "final_max_eaf": stats["max"],
        "final_mean_eaf": stats["mean"],
        "final_total_variants": stats["total"],
        "final_eaf_col_used": sample_column_dict["eaf_col"],
    })

    log_print("🎯 EAF harmonization complete.\n")

    # -------------------------------------------------------
    # WRITE LOG FILE
    # -------------------------------------------------------
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return df_checked, qc_info, sample_column_dict




'''
==============================================================
🧠 LOGIC SUMMARY: EAF / MAF HARMONIZATION PROCESS
==============================================================

1️⃣ Check if internal EAF column exists
    - If present → calculate proportion of variants with EAF > 0.6
    - If >25% → column is true EAF → compute zmaf = min(EAF, 1−EAF)
    - If ≤25% → treat as MAF → will be corrected using external EAF file

2️⃣ If EAF column missing or "NA"
    - Check if external EAF file (eaffile) exists
    - If not, use default reference file (default_eaf_file)

3️⃣ Load external EAF reference
    - Standardize columns: "chr", "pos", "a1", "a2" → renamed to match GWAS df
      (chr_col, pos_col, ea_col, oa_col)
    - Normalize alleles (uppercase), chromosomes, and types

4️⃣ Merge GWAS and EAF data
    - Merge once by allele orientation (A1/A2)
    - Merge again with flipped alleles (A2/A1), adjust EAF = 1−EAF
    - Concatenate both, then deduplicate by (chr, pos, A1, A2)

5️⃣ Compute minor allele frequency (MAF)
    - For each variant: zmaf = EAF if EAF ≤ 0.5 else 1−EAF

6️⃣ Record QC information
    - Variant counts before & after merging (direct, flipped, total)
    - File used, classification (EAF or MAF), and final variant count

7️⃣ Return
    - Harmonized dataframe with verified EAF and zmaf columns
    - QC dictionary summarizing all counts and processing steps
==============================================================

'''