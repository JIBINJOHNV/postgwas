import argparse
from pathlib import Path

from postgwas.utils.main import (
    validate_path,
    detect_total_memory_gb,
    validate_alphanumeric,
    validate_prefix_files,
)

# MAGMA covariate engine (step 6)
from postgwas.magmacovar.main import run_magma_covariates

# MAGMA gene-level analysis (step 5)
from postgwas.gene_assoc.magma_main import magma_analysis_pipeline


# ==========================================================
# DIRECT MODE ENGINE
# ==========================================================
def run_magma_covar_direct(args: argparse.Namespace) -> None:
    """
    Run ONLY the MAGMA covariate (gene-property) analysis step.
    """
    
    print("⚠ run_magma_covar_direct started")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Prefix = <outdir>/<sample_id>
    output_prefix = str(outdir / args.sample_id)

    if args.covar_model:
        print("⚠ WARNING: Using a custom covariate model (--covar_model).")
        print("  Ensure this matches your scientific intent.\n")

    # ---------------------------------------------------------
    # Call MAGMA covariates engine
    # ---------------------------------------------------------
    run_magma_covariates(
        magma_bin=str(args.magma),
        gene_results_file=str(args.magama_gene_assoc_raw),
        covariates_file=str(args.covariates),
        output_prefix=output_prefix,
        model=args.covar_model,
        direction=args.covar_direction,
        log_file=getattr(args, "log_file", None),
    )

    print("\n🎉 MAGMA Covariate Analysis (Direct Mode) Completed.")
    print(f"Results prefix: {output_prefix}\n")


# ==========================================================
# PIPELINE MODE ENGINE
# ==========================================================
def run_magma_covar_pipeline(args: argparse.Namespace) -> None:
    """
    Light MAGMA pipeline:
        1) MAGMA gene-level analysis
        2) MAGMA covariate analysis
    This pipeline assumes harmonisation, LD-blocking and formatting
    have been handled upstream or are not required.
    """

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Prefix for output files
    prefix_gene = str(outdir / f"{args.sample_id}_magma_gene")
    prefix_covar = str(outdir / f"{args.sample_id}_magma_covar")

    print("\n============================")
    print("🔬 Step 1: MAGMA Gene-Level Analysis")
    print("============================\n")

    # ---------------------------------------------------------
    # Run MAGMA gene-level analysis
    # ---------------------------------------------------------
    magma_gene_results = magma_analysis_pipeline(
        magma_bin=str(args.magma),
        plink_prefix=str(args.ref_panel),
        sumstats_vcf=str(args.vcf),
        gene_loc_file=str(args.gene_loc),
        out_prefix=prefix_gene,
        n_threads=args.nthreads,
        log_file=getattr(args, "log_file", None),
    )

    gene_raw_path = magma_gene_results.get("gene_raw")
    if not gene_raw_path:
        raise RuntimeError("MAGMA gene analysis did not produce a .genes.raw file.")

    print("✔ MAGMA gene analysis completed.")
    print(f"  Gene results: {gene_raw_path}\n")

    print("============================")
    print("🔬 Step 2: MAGMA Covariate (Gene-Property) Analysis")
    print("============================\n")

    # ---------------------------------------------------------
    # Run MAGMA covariate analysis
    # ---------------------------------------------------------
    run_magma_covariates(
        magma_bin=str(args.magma),
        gene_results_file=str(args.magama_gene_assoc_raw),
        covariates_file=str(args.covariates),
        output_prefix=prefix_covar,
        model=args.covar_model,
        direction=args.covar_direction,
        log_file=getattr(args, "log_file", None),
    )

    print("\n🎉 MAGMA Covariate Pipeline Completed Successfully")
    print(f"Gene-level results : {gene_raw_path}")
    print(f"Covariate output   : {prefix_covar}\n")
