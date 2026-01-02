import polars as pl
import shutil
import os
import logging
import multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import time

# Internal Imports
from postgwas.finemap.finemap.input_gen import (
    write_finemap_z_file, 
    write_snp_file,
    create_ldstore_master, 
    create_finemap_master
)
from postgwas.finemap.finemap.runner import (
    run_plink_extraction, 
    run_plink_to_bgen,
    run_bgen_indexing, 
    run_ldstore, 
    run_finemap_binary
)
from postgwas.finemap.finemap.parser import (
    parse_finemap_cred_file, 
    parse_finemap_config_file
)

# Configure Logger
logger = logging.getLogger("postgwas.finemap")

# ==============================================================
# 1. DEPENDENCY CHECK
# ==============================================================
def check_dependencies():
    """Ensures all required external binaries are available."""
    required_tools = ["plink2", "bgenix", "ldstore", "finemap"]
    missing = []
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    
    if missing:
        raise EnvironmentError(
            f"❌ Missing required tools: {', '.join(missing)}. "
            "Please install them and ensure they are in your PATH."
        )

# ==============================================================
# 2. WORKER FUNCTION (EXECUTES ON ONE LOCUS)
# ==============================================================
def process_single_locus(task):
    """
    Worker function to run the full FINEMAP pipeline for a single locus.
    """
    locus_id = task["locus_id"]
    out_dir = task["out_dir"]
    z_file_src = task["z_file"]     # Path to pre-sliced .z file
    snp_file_src = task["snp_file"] # Path to pre-sliced .snp file
    ld_ref_prefix = task["ld_ref"]
    n_samples = task["n_samples"]
    config = task["finemap_config"] # Dict of CLI args
    
    try:
        # Create unique directory for this locus
        locus_dir = out_dir / locus_id
        locus_dir.mkdir(parents=True, exist_ok=True)

        # --------------------------------------------------
        # A. Setup File Paths
        # --------------------------------------------------
        # Inputs (Copy to locus dir)
        z_file      = locus_dir / f"{locus_id}.z"
        snp_file    = locus_dir / f"{locus_id}.snp"
        
        # Intermediate (PLINK/BGEN)
        plink_out   = locus_dir / f"{locus_id}_plink"
        bgen_prefix = locus_dir / f"{locus_id}_genotypes"
        bgen_file   = locus_dir / f"{locus_id}_genotypes.bgen"
        bgi_file    = locus_dir / f"{locus_id}_genotypes.bgen.bgi"
        
        # Intermediate (LDstore)
        ldstore_master = locus_dir / f"{locus_id}.ldstore.master"
        bcor_file      = locus_dir / f"{locus_id}.bcor"
        ld_matrix      = locus_dir / f"{locus_id}.ld"
        
        # Outputs (FINEMAP)
        master_file = locus_dir / f"{locus_id}.master"
        config_file = locus_dir / f"{locus_id}.config" # Output by finemap
        cred_file   = locus_dir / f"{locus_id}.cred"   # Output by finemap
        log_file    = locus_dir / f"{locus_id}.log"    # Output by finemap

        # Copy pre-sliced inputs to working dir
        shutil.copy(z_file_src, z_file)
        shutil.copy(snp_file_src, snp_file)

        # --------------------------------------------------
        # B. LD Generation Pipeline
        # --------------------------------------------------
        
        # 1. PLINK Extraction (Ref -> BED)
        try:
            run_plink_extraction(ld_ref_prefix, snp_file, plink_out)
        except RuntimeError as e:
            if "Exit status 13" in str(e) or "No variants" in str(e):
                return f"⚠️ {locus_id}: Skipped (0 variants matched Reference Panel IDs)"
            raise e

        # 2. PLINK Export (BED -> BGEN 8-bit)
        run_plink_to_bgen(plink_out, bgen_prefix)

        # 3. Index BGEN
        run_bgen_indexing(bgen_file)

        # 4. Determine LD Reference Sample Size
        fam_file = Path(f"{ld_ref_prefix}.fam")
        if not fam_file.exists():
            raise FileNotFoundError(f"Reference .fam file not found at: {fam_file}\nCheck your --finemap_ld_ref path.")
            
        with open(fam_file, 'rb') as f:
            n_ld_samples = sum(1 for _ in f)

        # 5. Create LDstore Master File
        create_ldstore_master(
            master_path=ldstore_master,
            z_file=z_file,
            bgen_file=bgen_file,
            bgi_file=bgi_file,
            bcor_file=bcor_file,
            ld_matrix=ld_matrix,
            n_samples=n_ld_samples
        )

        # 6. Run LDstore
        run_ldstore(ldstore_master, n_threads=1) 

        # --------------------------------------------------
        # C. Run FINEMAP (With Logic Fix)
        # --------------------------------------------------
        
        # FIX: Check actual number of SNPs available
        with open(snp_file, 'r') as f:
            # Count non-empty lines
            n_snps_available = sum(1 for line in f if line.strip())
            
        # Create a local config copy so we don't mess up other workers
        local_config = config.copy()
        
        # Dynamic Adjustment: n_causal_snps cannot exceed total SNPs
        if n_snps_available < local_config["n_causal_snps"]:
            local_config["n_causal_snps"] = n_snps_available
            
        # 1. Create FINEMAP Master File
        create_finemap_master(
            master_path=master_file,
            z_file=z_file,
            ld_file=ld_matrix,
            snp_file=snp_file,
            config_file=config_file,
            cred_file=cred_file,
            log_file=log_file,
            n_samples=n_samples
        )

        # 2. Execute Binary using the ADJUSTED config
        success, msg = run_finemap_binary(master_file, local_config, n_threads=1)

        if not success:
            return f"⚠️ {locus_id}: {msg}"

        # --------------------------------------------------
        # D. Cleanup
        # --------------------------------------------------
        for f in [plink_out.with_suffix(".bed"), plink_out.with_suffix(".bim"), 
                  plink_out.with_suffix(".fam"), bgen_file, bgi_file, bcor_file]:
            if f.exists():
                f.unlink()

        return {
            "status": "success",
            "locus_id": locus_id,
            "cred_file": str(cred_file),
            "config_file": str(config_file)
        }

    except Exception as e:
        return f"❌ {locus_id} Failed: {str(e)}"


# ==============================================================
# 3. MAIN DRIVER (CALLED BY CLI)
# ==============================================================
def run_finemap_pipeline(args):
    """
    Main driver for the FINEMAP pipeline.
    """
    start_time = time()
    check_dependencies()
    
    # 1. Setup Paths
    outdir = Path(args.outdir).resolve()
    temp_dir = outdir / "temp_inputs"
    results_dir = outdir / "loci_data"
    final_output_dir = outdir / "final_results"
    
    for d in [outdir, temp_dir, results_dir, final_output_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # 2. Load Data
    logger.info(f"Loading Locus File: {args.locus_file}")
    loci_df = pl.read_csv(args.locus_file, separator="\t")
    print(f"total loci is {loci_df.shape}")
    logger.info(f"Loading Summary Statistics: {args.finemap_in_files}")
    sumstats_df = pl.read_csv(args.finemap_in_files, separator="\t")
    # 3. MHC Filter
    if args.finemap_skip_mhc and not args.finemap_include_mhc:
        logger.info(f"Filtering MHC Region (Chr{args.finemap_mhc_chrom})")
        loci_df = loci_df.filter(
            ~((pl.col("CHR") == args.finemap_mhc_chrom) & 
              (pl.col("START") < args.finemap_mhc_end) & 
              (pl.col("END") > args.finemap_mhc_start))
        )
    logger.info("Harmonizing Variant IDs to Reference Format (CHR_POS_A2_A1)...")
    # Ensure chromosome is string and clean (remove 'chr' if present)
    sumstats_df = sumstats_df.with_columns(
        pl.col("chromosome").cast(pl.Utf8).str.replace("chr", "")
    )
    # =========================================================
    # 5. INTERSECTION WITH REFERENCE PANEL (New Block)
    # =========================================================
    # PLINK Prefix Logic
    raw_ref = str(args.finemap_ld_ref)
    if raw_ref.endswith(".bed") or raw_ref.endswith(".bim") or raw_ref.endswith(".fam"):
         ld_ref_prefix = str(Path(raw_ref).with_suffix(""))
    else:
         ld_ref_prefix = raw_ref
    finemap_ld_ref_bim = f"{ld_ref_prefix}.bim"
    if not Path(finemap_ld_ref_bim).exists():
        raise FileNotFoundError(f"Reference BIM file not found at: {finemap_ld_ref_bim}")
    logger.info(f"Reading LD reference BIM: {finemap_ld_ref_bim}")
    ld_ref_df = pl.read_csv(
        finemap_ld_ref_bim,
        separator="\t",
        has_header=False,
        columns=[1], # BIM col 2 is variant ID
        new_columns=["variant_id"],
    )
    # AF → MAF transformation (if not already MAF)
    sumstats_df = sumstats_df.with_columns(
        pl.when(pl.col("maf") <= 0.5)
        .then(pl.col("maf"))
        .otherwise(1.0 - pl.col("maf"))
        .alias("maf")
    )
    # Deterministic order
    sumstats_df = sumstats_df.sort(
        by=["chromosome", "position", "allele1", "allele2"]
    )
    # Filter to LD reference intersection
    logger.info("Intersecting with Reference Panel...")
    sumstats_filt = sumstats_df.join(
        ld_ref_df,
        left_on="rsid",
        right_on="variant_id",
        how="inner", 
    )
    logger.info(f"Variants remaining after intersection: {sumstats_filt.height} (from {sumstats_df.height})")
    if sumstats_filt.is_empty():
        raise ValueError("❌ No variants remain after intersection! Check your ID formats match the BIM file.")
    # =========================================================
    logger.info("Preparing locus-specific input files...")
    tasks = []
    finemap_config = {
        "n_causal_snps": args.n_causal_snps,
        "prob_cred_set": args.prob_cred_set,
        "n_iter": args.n_iter,
        "n_conv_sss": args.n_conv_sss,
        "prob_conv_sss_tol": args.prob_conv_sss_tol,
        "prior_k": args.prior_k,
        "std_effects": args.std_effects,
    }
    # Loop through loci
    chrom_col = "CHR" if "CHR" in loci_df.columns else "chromosome"
    print(f"total loci is {loci_df.shape}")
    for row in loci_df.iter_rows(named=True):
        chrom = str(row[chrom_col]).replace("chr", "")
        start = row["START"]
        end = row["END"]
        locus_id = f"chr{chrom}_{start}_{end}"
        # Slice Sumstats (Using FILTERED data now)
        locus_ss = sumstats_filt.filter(
            (pl.col("chromosome") == chrom) &
            (pl.col("position") >= start) &
            (pl.col("position") <= end)
        )
        if locus_ss.height < 1:
            # logger.warning(f"Skipping {locus_id}: Too few variants ({locus_ss.height})")
            continue
        # Get Sample Size (NEF)
        if "NEF" in locus_ss.columns:
            n_samples = int(locus_ss.select(pl.col("NEF").drop_nulls().median()).item())
        else:
            logger.warning(f"{locus_id}: NEF missing, using default 10,000")
            n_samples = 10000
        # Create Temp Files (Main Thread)
        z_outfile = temp_dir / f"{locus_id}.z"
        snp_outfile = temp_dir / f"{locus_id}.snp"
        write_finemap_z_file(locus_ss, z_outfile)
        write_snp_file(locus_ss, snp_outfile)
        tasks.append({
            "locus_id": locus_id,
            "z_file": z_outfile,
            "snp_file": snp_outfile,
            "out_dir": results_dir,
            "ld_ref": ld_ref_prefix, 
            "n_samples": n_samples,
            "finemap_config": finemap_config
        })
    logger.info(f"Prepared {len(tasks)} loci for parallel analysis.")    
    
    # =========================================================
    # 5. EXECUTE PARALLEL (MEMORY AWARE)
    # =========================================================
    # Calculate Memory Limit
    MIN_RAM_PER_WORKER_GB = 14
    
    # Get Total RAM in GB (Cross-platform safe)
    try:
        total_ram_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        total_ram_gb = total_ram_bytes / (1024.**3)
    except (ValueError, AttributeError):
        # Fallback if sysconf fails (rare on Linux/Mac)
        total_ram_gb = 16 # Assume 16GB as conservative default
    
    max_workers_ram = int(total_ram_gb // MIN_RAM_PER_WORKER_GB)
    max_workers_cpu = getattr(args, "threads", max(1, os.cpu_count() - 2))
    
    # Final workers = Min(CPU limit, RAM limit), but at least 1
    max_workers = max(1, min(max_workers_cpu, max_workers_ram))
    
    
    successful_results = []
    print(f"Total tasks: {len(tasks)}")
    # Use 'spawn' for stability with subprocesses
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as executor:
        futures = [executor.submit(process_single_locus, t) for t in tasks]
        for i, fut in enumerate(as_completed(futures)):
            result = fut.result()
            # Handle Return
            if isinstance(result, dict) and result.get("status") == "success":
                successful_results.append(result)
                if (i + 1) % 5 == 0:
                    logger.info(f"[{i+1}/{len(tasks)}] Locus {result['locus_id']} finished.")
            else:
                # Only log real errors, not skips
                level = logging.WARNING if "Skipped" in str(result) else logging.ERROR
                logger.log(level, f"[{i+1}/{len(tasks)}] {result}")
    # # 6. Aggregation / Post-Processing
    # logger.info("Aggregating results...")
    # all_cred = []
    # all_config = []
    # for res in successful_results:
    #     locus_id = res['locus_id']
    #     cred_path = Path(res['cred_file'])
    #     config_path = Path(res['config_file'])
    #     # Use parsers from parser.py
    #     df_cred = parse_finemap_cred_file(cred_path, locus_id)
    #     if not df_cred.is_empty():
    #         all_cred.append(df_cred)
    #     df_config = parse_finemap_config_file(config_path, locus_id)
    #     if not df_config.is_empty():
    #         all_config.append(df_config)
    # # 7. Write Final Output
    # if all_cred:
    #     final_cred_df = pl.concat(all_cred)
    #     out_cred = final_output_dir / "finemap_results_cred.tsv"
    #     final_cred_df.write_csv(out_cred, separator="\t")
    #     logger.info(f"Saved SNP-level results: {out_cred}")
    # else:
    #     logger.warning("No SNP-level results found!")
    # if all_config:
    #     final_config_df = pl.concat(all_config)
    #     out_config = final_output_dir / "finemap_results_config.tsv"
    #     final_config_df.write_csv(out_config, separator="\t")
    #     logger.info(f"Saved Config results: {out_config}")
    # # 8. Cleanup Temp Inputs
    # if not getattr(args, "keep_temp", False):
    #     shutil.rmtree(temp_dir, ignore_errors=True)
    # logger.info(f"✅ Pipeline completed in {time() - start_time:.2f} seconds.")