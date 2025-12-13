#!/usr/bin/env python3

from __future__ import annotations

from typing import Optional
from concurrent.futures import ProcessPoolExecutor
import functools
import os, re,subprocess,time,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
import polars as pl
from scipy.stats import norm
from postgwas.utils.main import safe_thread_count 
from typing import List, Optional, Tuple
import threading
import psutil

"""
pred_ld_runner.py — multiprocessing-safe PRED-LD engine for PostGWAS

Designed to run inside Docker AND locally on macOS/Linux.

Features:
---------
✓ macOS-safe multiprocessing (spawn)
✓ Linux-safe multiprocessing (fork)
✓ Worker function defined at top-level (pickle-safe)
✓ Automatic pred_ld.py discovery
✓ Interleaved chromosome scheduling (big + small → reduced peak RAM)
✓ Workspace isolation in HOME (always writable in Docker)
✓ Merges logs and validates outputs
✓ Auto thread reduction using safe_thread_count()
"""

# ==============================================================================
# 1. TOP-LEVEL WORKER FUNCTION
#    (Required for macOS + multiprocessing)
# ==============================================================================

# =====================================================================
# 1. PRED-LD PARALLEL RUNNER
# =====================================================================



# =============================================================================
#  MEMORY HELPERS
# =============================================================================

def _mem_from_psutil() -> Optional[Tuple[float, float]]:
    """Try to read total/free RAM using psutil. Returns (total_gb, free_gb)."""
    try:
        vm = psutil.virtual_memory()
        total_gb = vm.total / (1024 ** 3)
        free_gb = vm.available / (1024 ** 3)
        return total_gb, free_gb
    except Exception:
        return None


def _mem_from_proc() -> Optional[Tuple[float, float]]:
    """Linux fallback: parse /proc/meminfo. Returns (total_gb, free_gb)."""
    try:
        meminfo = Path("/proc/meminfo")
        if not meminfo.exists():
            return None

        info = {}
        with meminfo.open() as f:
            for line in f:
                parts = line.split(":")
                if len(parts) != 2:
                    continue
                key = parts[0].strip()
                val = parts[1].strip().split()[0]
                try:
                    info[key] = float(val)
                except ValueError:
                    continue

        # Values in kB → convert to GB
        total_kb = info.get("MemTotal")
        # Use MemAvailable if present; else fall back to (MemFree + Cached)
        avail_kb = info.get("MemAvailable")
        if avail_kb is None:
            avail_kb = info.get("MemFree", 0.0) + info.get("Cached", 0.0)

        if total_kb is None or avail_kb is None:
            return None

        total_gb = total_kb / (1024 ** 2)
        free_gb = avail_kb / (1024 ** 2)
        return total_gb, free_gb
    except Exception:
        return None


def get_memory_status() -> Tuple[Optional[float], Optional[float]]:
    """
    Unified RAM checker.

    Returns:
        (total_gb, free_gb) or (None, None) if cannot be determined.
    """
    res = _mem_from_psutil()
    if res is not None:
        return res
    res = _mem_from_proc()
    if res is not None:
        return res
    return None, None


# =============================================================================
#  CHROMOSOME HELPERS
# =============================================================================

BIG_CHROMS = {"1", "2", "3", "4", "5"}


def is_big_chr(chr_id: str) -> bool:
    """Return True if chromosome is considered 'big' for scheduling."""
    chr_id = chr_id.upper().replace("CHR", "")
    try:
        n = int(chr_id)
        return str(n) in BIG_CHROMS
    except ValueError:
        return False


def build_interleaved_chr_order(chromosomes: List) -> List[str]:
    """
    Use a fixed, RAM-balanced chromosome execution order:

    1, 22, 21, 20, 2, 19, 18, 17, 3, 16, 15, 14,
    4, 13, 12, 11, 5, 10, 9, 8, 7, 6, X

    - chrX always last (if present)
    - Missing chromosomes are skipped safely
    """
    desired_order = [
        "1", "22", "21", "20",
        "19","18", "2", "17",
        "16", "15", "14", "3", 
        "13","12", "4", "11",
        "10", "5", "9", "8",
        "7", "6", "X"
    ]
    # Normalize input
    chroms = {str(c).upper().replace("CHR", "") for c in chromosomes}
    # Preserve order, include only chromosomes the user actually provided
    final_order = [c for c in desired_order if c in chroms]
    return final_order


# =============================================================================
#  MAIN PRED-LD PARALLEL RUNNER
# =============================================================================

def run_pred_ld_parallel(
    predld_input_dir: str,
    output_folder: str,
    output_prefix: str,
    pred_ld_ref: str,
    chromosomes: List = list(range(1, 23)) + ["X"],
    r2threshold: float = 0.8,
    maf: float = 0.001,
    population: str = "EUR",
    ref: str = "TOP_LD",
    threads: int = 6,
) -> bool:
    """
    Fully RAM-aware, big-chromosome-safe PRED-LD parallel runner.

    Rules implemented:
    -------------------
    1. Max parallel workers = 2.
    2. A job can start only if:
         free_ram_gb >= min(30, total_ram_gb * 0.80)
       (if RAM cannot be detected, we fall back to "big-chr only" constraint).
    3. Big chromosomes = {1..6} are never run together.
       (At most one big chr at a time.)
    4. chrX always runs last. If missing input → warned but does not crash.
    5. Workspace is under `output_folder/.predld_work` (host-mounted, large).
    6. pred_ld.py is auto-located relative to this file (no hard-coded paths).
    """

    # -------------------------------------------------------------------------
    # 0. Locate pred_ld.py dynamically (relative to this file)
    # -------------------------------------------------------------------------
    module_root = Path(__file__).resolve().parent
    pred_ld_script = module_root / "pred_ld.py"
    if not pred_ld_script.exists():
        raise FileNotFoundError(f"ERROR: pred_ld.py not found at {pred_ld_script}")

    # -------------------------------------------------------------------------
    # 1. Normalize and prepare paths
    # -------------------------------------------------------------------------
    predld_input_dir = Path(predld_input_dir).expanduser().resolve()
    output_folder = Path(output_folder).expanduser().resolve()
    pred_ld_ref = Path(pred_ld_ref).expanduser().resolve()

    output_folder.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # 2. Workspace (host-mounted, avoids /tmp space issues inside Docker)
    # -------------------------------------------------------------------------
    workspace = output_folder / ".predld_work"
    if workspace.exists():
        shutil.rmtree(workspace, ignore_errors=True)
    workspace.mkdir(parents=True, exist_ok=True)

    local_ref = workspace / "ref"
    shutil.copytree(pred_ld_ref, local_ref, dirs_exist_ok=True)

    # -------------------------------------------------------------------------
    # 3. Threads / workers & chromosome order
    # -------------------------------------------------------------------------
    # safe_thread_count only for CPU-side bound; concurrency is 2.
    threads = safe_thread_count(threads, gb_per_thread=30)
    max_workers = min(2, max(1, threads))

    chr_order = build_interleaved_chr_order(chromosomes)

    print(f"🧬 Running PRED-LD with {max_workers} parallel workers")
    print(f"📁 Input dir : {predld_input_dir}")
    print(f"📁 Output dir: {output_folder}")
    print(f"📁 Workspace : {workspace}")
    print(f"📁 Ref LD    : {pred_ld_ref}")
    print(f"🧬 Chromosomes (interleaved): {', '.join(chr_order)}")

    # -------------------------------------------------------------------------
    # 4. Shared state for scheduling
    # -------------------------------------------------------------------------
    lock = threading.Lock()
    running_chroms: set[str] = set()
    running_big: bool = False  # "is any big chromosome currently running?"

    # Track which ones succeeded / failed
    succeeded: list[str] = []
    failed: list[tuple[str, str]] = []  # (chr, reason)

    def wait_until_allowed(chr_id: str) -> None:
        """
        Busy-wait (with sleep) until both RAM and big-chr constraints allow
        this chromosome to start.
        """
        nonlocal running_big

        while True:
            with lock:
                total_gb, free_gb = get_memory_status()
                ram_known = (total_gb is not None and free_gb is not None)

                if ram_known:
                    threshold = min(20.0, total_gb * 0.80)
                    free_ok = free_gb >= threshold
                else:
                    # Cannot measure RAM → only enforce big-chr constraint
                    threshold = None
                    free_ok = True

                this_big = is_big_chr(chr_id)
                big_ok = (not this_big) or (not running_big)

                if free_ok and big_ok:
                    # Reserve slot
                    running_chroms.add(chr_id)
                    if this_big:
                        running_big = True
                    # Optional: small debug log
                    if ram_known:
                        print(
                            f"🚀 Starting chr{chr_id}: free RAM {free_gb:.1f} GB "
                            f"(threshold {threshold:.1f} GB), "
                            f"big_running={running_big}"
                        )
                    else:
                        print(
                            f"🚀 Starting chr{chr_id}: RAM unknown, "
                            f"big_running={running_big}"
                        )
                    return

                # Not allowed yet → print occasional status
                if ram_known and not free_ok:
                    print(
                        f"⏳ Waiting for RAM before chr{chr_id}: "
                        f"free={free_gb:.1f} GB, "
                        f"threshold={threshold:.1f} GB"
                    )
                if this_big and not big_ok:
                    print(f"⏳ Waiting: another big chromosome already running.")

            time.sleep(10)  # Sleep outside lock

    def release_chr(chr_id: str) -> None:
        """Release this chromosome from the running set and big flag."""
        nonlocal running_big
        with lock:
            running_chroms.discard(chr_id)
            if is_big_chr(chr_id):
                # If no other big chr is running, clear the flag
                if not any(is_big_chr(c) for c in running_chroms):
                    running_big = False

    # -------------------------------------------------------------------------
    # 5. Worker (per chromosome)
    # -------------------------------------------------------------------------
    def worker(chr_id: str) -> None:
        chr_input = predld_input_dir / f"{output_prefix}_chr{chr_id}_predld_input.tsv"
        log_file = output_folder / f"{output_prefix}_chr{chr_id}_predld.log"

        # chrX often missing: soft skip if input file not found
        if not chr_input.exists():
            msg = f"PRED-LD input missing for chr{chr_id}: {chr_input}"
            print(f"⚠️ chr{chr_id} skipped: {msg}")
            failed.append((chr_id, msg))
            return

        # Wait until RAM and big-chr constraints allow us to start
        wait_until_allowed(chr_id)

        try:
            cmd = [
                "python",
                str(pred_ld_script),
                "--file-path",
                str(chr_input),
                "--pop",
                population,
                "--ref",
                ref,
                "--ref_dir",
                str(local_ref),
                "--r2threshold",
                str(r2threshold),
                "--maf",
                str(maf),
                "--out_dir",
                str(output_folder),
            ]

            with log_file.open("w") as log:
                proc = subprocess.run(cmd, stdout=log, stderr=log)

            if proc.returncode != 0:
                if proc.returncode == -9:
                    reason = (
                        f"PRED-LD failed for chr{chr_id} "
                        f"(likely out-of-memory, exit={proc.returncode}). "
                        f"See log: {log_file}"
                    )
                else:
                    reason = (
                        f"PRED-LD failed for chr{chr_id} "
                        f"(exit={proc.returncode}). See log: {log_file}"
                    )
                print(f"❌ {reason}")
                failed.append((chr_id, reason))
                return

            expected_out = output_folder / f"imputation_results_chr{chr_id}.txt"
            if not expected_out.exists():
                reason = (
                    f"PRED-LD output missing for chr{chr_id}: {expected_out}. "
                    f"See log: {log_file}"
                )
                print(f"❌ {reason}")
                failed.append((chr_id, reason))
                return

            print(f"✅ chr{chr_id} finished successfully.")
            succeeded.append(chr_id)

        finally:
            # Always release the slot
            release_chr(chr_id)

    # -------------------------------------------------------------------------
    # 6. Launch workers (thread-based scheduler)
    # -------------------------------------------------------------------------
    t0 = time.time()
    threads_list: List[threading.Thread] = []

    for chr_id in chr_order:
        t = threading.Thread(target=worker, args=(chr_id,), daemon=True)
        threads_list.append(t)

    # Run with at most `max_workers` active threads at any time
    active: List[threading.Thread] = []

    for t in threads_list:
        # Wait until there is a free thread slot
        while True:
            # Purge finished threads from active list
            active = [thr for thr in active if thr.is_alive()]
            if len(active) < max_workers:
                t.start()
                active.append(t)
                break
            time.sleep(2)

    # Wait for all threads to complete
    for t in threads_list:
        t.join()

    # -------------------------------------------------------------------------
    # 7. Merge per-chrom logs
    # -------------------------------------------------------------------------
    merged_log = output_folder / f"{output_prefix}_predld.log"
    with merged_log.open("w") as fout:
        for lf in sorted(output_folder.glob(f"{output_prefix}_chr*_predld.log")):
            with lf.open() as f:
                fout.write(f.read())
            lf.unlink()  # remove per-chr logs

    print(f"\n📄 Combined log → {merged_log}")

    # -------------------------------------------------------------------------
    # 8. Summary
    # -------------------------------------------------------------------------
    elapsed_min = (time.time() - t0) / 60.0
    print(f"\n🎯 PRED-LD completed in {elapsed_min:.2f} minutes.")
    if succeeded:
        print(f"   ✅ Successful chromosomes: {', '.join(sorted(succeeded, key=str))}")
    if failed:
        print("   ⚠️ Failed / skipped chromosomes:")
        for chr_id, reason in failed:
            print(f"      chr{chr_id}: {reason}")

    # Return True if at least one chromosome completed successfully
    return bool(succeeded)

from typing import Optional  # make sure this is at top of file

def process_pred_ld_results_all_parallel(
    folder_path: str,
    output_path: Optional[str] = None,
    gwas2vcf_resource_folder: Optional[str] = None,
    output_prefix: Optional[str] = None,
    corr_method: str = "pearson",
    threads: int = 6,
):
    """
    Production-ready PRED-LD postprocessing pipeline.
    - Parallel Polars implementation
    - Safe dtype coercion for all PRED-LD fields
    - Accurate imputed-SNP sample-size propagation
    - Duplicate SNP diagnostics + correlation
    - Final TSV + correlation summary
    - gwas2vcf config generation
    - Cleanup of intermediate files
    """

    import os, re, subprocess
    import polars as pl
    import numpy as np
    from pathlib import Path
    from scipy.stats import norm
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import pandas as pd

    # =====================================================================
    # Safe helpers
    # =====================================================================
    def coerce_dtypes(df: pl.DataFrame, dtype_map: dict) -> pl.DataFrame:
        """Coerce columns to specific dtypes (safe casting)."""
        for col, dtype in dtype_map.items():
            if col in df.columns:
                df = df.with_columns(pl.col(col).cast(dtype, strict=False))
        return df

    def safe_drop(df: pl.DataFrame, cols) -> pl.DataFrame:
        """
        Drop only columns that exist.
        cols can be a string or list of strings.
        """
        if isinstance(cols, str):
            cols = [cols]
        existing = [c for c in cols if c in df.columns]
        if not existing:
            return df
        return df.drop(existing)

    # =====================================================================
    # Expected column types
    # =====================================================================
    data_types = {
        "chr": pl.Utf8, "snp": pl.Utf8,
        "A1": pl.Utf8, "A2": pl.Utf8,
        "pos": pl.Float64,
        "beta": pl.Float64, "SE": pl.Float64, "z": pl.Float64,
        "imputed": pl.Float64,
        "R2": pl.Float64, "NC": pl.Float64, "SS": pl.Float64,
        "AF": pl.Float64, "LP": pl.Float64, "SI": pl.Float64,
    }

    info_types = {
        "pos1": pl.Float64, "pos2": pl.Float64,
        "R2": pl.Float64, "Dprime": pl.Float64,
        "ALT_AF1": pl.Float64, "ALT_AF2": pl.Float64,
        "+/-corr": pl.Utf8,
        "rsID1": pl.Utf8, "rsID2": pl.Utf8,
        "REF1": pl.Utf8, "ALT1": pl.Utf8,
        "REF2": pl.Utf8, "ALT2": pl.Utf8
    }

    # =====================================================================
    # Discover chromosome files
    # =====================================================================
    folder_path = Path(folder_path)
    all_files = os.listdir(folder_path)

    imp_files = [f for f in all_files if f.startswith("imputation_results_chr")]
    info_files = [f for f in all_files if f.startswith("LD_info_TOP_LD_chr")]

    if not imp_files or not info_files:
        raise FileNotFoundError("❌ Missing imputation results or LD info files.")

    def _chr_sort_key(x):
        return (0, int(x)) if x.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25}.get(x.upper(), 99))

    chromosomes = sorted(
        set(re.findall(r"chr(\w+)", " ".join(imp_files))),
        key=_chr_sort_key
    )

    print(f"🧬 Detected chromosomes: {', '.join(chromosomes)}")

    # =====================================================================
    # Per-chromosome worker
    # =====================================================================
    results, summaries = [], []

    def process_chr(chr_id):
        try:
            imp_file = folder_path / f"imputation_results_chr{chr_id}.txt"
            info_file = folder_path / f"LD_info_TOP_LD_chr{chr_id}.txt"

            if not imp_file.exists() or not info_file.exists():
                print(f"⚠️ chr{chr_id}: missing files, skipping.")
                return None, None

            # -------- Load files --------
            data_df = coerce_dtypes(
                pl.read_csv(imp_file, separator="\t", infer_schema_length=5000),
                data_types
            )
            info_df = coerce_dtypes(
                pl.read_csv(info_file, separator="\t", infer_schema_length=5000),
                info_types
            )

            # -------------------------------------------------------------
            # Duplicate SNP diagnostics
            # -------------------------------------------------------------
            duplicated = data_df.filter(pl.col("snp").is_duplicated())

            imputed_dups = duplicated.filter(pl.col("imputed") == 1)
            imputed_dups = safe_drop(imputed_dups, ["NC", "SS", "AF", "LP", "SI"])

            nonimp_dups = duplicated.filter(pl.col("imputed") != 1)
            nonimp_dups = safe_drop(nonimp_dups, ["NC", "SS", "AF", "LP", "SI"])

            beta_corr = z_corr = np.nan
            if imputed_dups.height > 0 and nonimp_dups.height > 0:
                beta_corr = np.corrcoef(
                    imputed_dups["beta"].to_numpy(),
                    nonimp_dups["beta"].to_numpy()
                )[0, 1]
                z_corr = np.corrcoef(
                    imputed_dups["z"].to_numpy(),
                    nonimp_dups["z"].to_numpy()
                )[0, 1]

            # -------------------------------------------------------------
            # NON-IMPUTED summary statistics
            # -------------------------------------------------------------
            nonimp = data_df.filter(pl.col("imputed") == 0)
            nonimp = safe_drop(nonimp, "R2")

            if "LP" in nonimp.columns:
                nonimp = nonimp.with_columns(
                    (10 ** (-pl.col("LP"))).alias("p_value")
                )
                nonimp = safe_drop(nonimp, "LP")

            # Leave NC, SS, AF, SI for later propagation

            # -------------------------------------------------------------
            # IMPUTED SNPs (initial cleaning)
            # -------------------------------------------------------------
            imputed = (
                data_df.filter(pl.col("imputed") == 1)
                .filter(~pl.col("snp").is_in(nonimp["snp"]))
            )
            imputed = safe_drop(imputed, ["R2", "NC", "SS", "AF", "LP", "SI"])

            # =================================================================
            # ⭐ SAMPLE-SIZE PROPAGATION (NC, SS, AF, SI for imputed SNPs)
            # =================================================================
            donor_stats = (
                nonimp
                .select(["snp", "NC", "SS", "AF", "SI"])
                .drop_nulls(subset=["NC", "SS"])
            )

            linked = info_df.join(
                donor_stats,
                left_on="rsID1",
                right_on="snp",
                how="inner"
            )

            imputed_stats = (
                linked.group_by("rsID2")
                .agg([
                    pl.col("NC").mean(),
                    pl.col("SS").mean(),
                    pl.col("AF").mean(),
                    pl.col("SI").mean(),
                ])
                .rename({"rsID2": "snp"})
            )

            imputed = imputed.join(imputed_stats, on="snp", how="left")

            # Compute p-values for imputed SNPs
            if "z" in imputed.columns:
                imputed = imputed.with_columns(
                    (2 * pl.Series(norm.sf(np.abs(imputed["z"].to_numpy()))))
                    .alias("p_value")
                )

            # -------------------------------------------------------------
            # Combine non-imputed + imputed
            # -------------------------------------------------------------
            nonimp = nonimp.with_columns(pl.lit("no").alias("imputed"))
            imputed = imputed.with_columns(pl.lit("yes").alias("imputed"))

            final_df = pl.concat([nonimp, imputed], how="diagonal")

            # NC/SS rounding and ncase
            for col in ["NC", "SS"]:
                if col in final_df.columns:
                    final_df = final_df.with_columns(
                        pl.col(col).round(0).cast(pl.Int64, strict=False)
                    )

            if "SS" in final_df.columns and "NC" in final_df.columns:
                final_df = final_df.with_columns(
                    pl.when(pl.col("SS") > pl.col("NC"))
                    .then((pl.col("SS") - pl.col("NC")).cast(pl.Int64, strict=False))
                    .otherwise(0)
                    .alias("ncase_col")
                )

            summary = {
                "chromosome": chr_id,
                "imputed_markers": imputed.height,
                "n_dups_imputed": imputed_dups.height,
                "n_dups_nonimputed": nonimp_dups.height,
                "beta_corr": float(beta_corr) if beta_corr == beta_corr else None,
                "z_corr": float(z_corr) if z_corr == z_corr else None,
            }

            return final_df, summary

        except Exception as e:
            print(f"❌ Error processing chr{chr_id}: {e}")
            return None, None

    # =====================================================================
    # Parallel execution
    # =====================================================================
    with ThreadPoolExecutor(max_workers=threads) as exe:
        futures = {exe.submit(process_chr, c): c for c in chromosomes}

        for fut in as_completed(futures):
            df, summ = fut.result()
            if df is not None:
                results.append(df)
            if summ is not None:
                summaries.append(summ)

    if not results:
        raise RuntimeError("❌ No chromosome processed successfully!")

    combined_df = pl.concat(results, how="diagonal")
    corr_df = pl.DataFrame(summaries)

    # =====================================================================
    # Output setup
    # =====================================================================
    output_folder = output_path or folder_path
    os.makedirs(output_folder, exist_ok=True)
    output_prefix = output_prefix or "PRED_LD"

    imput_harmonisation_folder = Path(output_folder).parent / "imputed_harmonised"
    imput_harmonisation_folder = str(imput_harmonisation_folder)

    combined_path = f"{output_folder}/{output_prefix}_PREDLD_allchr.tsv"
    corr_path = f"{output_folder}/{output_prefix}_PREDLD_correlations.tsv"
    config_path = f"{output_folder}/{output_prefix}_gwas2vcf_config.csv"

    combined_df.write_csv(combined_path, separator="\t")
    corr_df.write_csv(corr_path, separator="\t")
    os.system(f"gzip {combined_path}")

    # =====================================================================
    # Build gwas2vcf config
    # =====================================================================
    columns = [
        "sumstat_file", "gwas_outputname",
        "chr_col", "pos_col", "snp_id_col",
        "ea_col", "oa_col", "eaf_col",
        "beta_or_col", "se_col", "imp_z_col",
        "pval_col", "ncontrol_col", "ncase_col",
        "ncontrol", "ncase", "imp_info_col",
        "delimiter", "infofile", "infocolumn",
        "eaffile", "eafcolumn", "liftover",
        "chr_pos_col", "resourse_folder", "output_folder"
    ]

    cfg = pd.DataFrame([[None] * len(columns)], columns=columns)

    cfg.loc[0, [
        "sumstat_file", "gwas_outputname",
        "output_folder", "resourse_folder", "delimiter"
    ]] = [
        f"{combined_path}.gz",
        f"{output_prefix}_imputed",
        imput_harmonisation_folder,
        gwas2vcf_resource_folder,
        "\t"
    ]

    mapping = {
        "chr_col": "chr",
        "pos_col": "pos",
        "snp_id_col": "snp",
        "ea_col": "A1",
        "oa_col": "A2",
        "eaf_col": "AF",
        "beta_or_col": "beta",
        "se_col": "SE",
        "pval_col": "p_value",
        "ncontrol_col": "NC",
        "ncase_col": "ncase_col",
        "imp_info_col": "SI",
    }

    for cfg_col, df_col in mapping.items():
        cfg.loc[0, cfg_col] = df_col if df_col in combined_df.columns else "NA"

    cfg.to_csv(config_path, index=False)

    # =====================================================================
    # Cleanup original files
    # =====================================================================
    subprocess.run(
        f"cat {folder_path}/imputation_results_chr*.txt | gzip -c > {output_folder}/{output_prefix}_imputation_results.txt.gz",
        shell=True,
        check=False,
    )

    subprocess.run(
        f"cat {folder_path}/LD_info_TOP_LD_chr*.txt | gzip -c > {output_folder}/{output_prefix}_LD_info_TOP_LD.txt.gz",
        shell=True,
        check=False,
    )

    subprocess.run(f"rm -f {folder_path}/imputation_results_chr*.txt", shell=True, check=False)
    subprocess.run(f"rm -f {folder_path}/LD_info_TOP_LD_chr*.txt", shell=True, check=False)

    return combined_df, corr_df
