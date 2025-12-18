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

# ==========================================================
# DIRECT MODE ENGINE
# ==========================================================
def run_magma_covar_direct(args: argparse.Namespace,ctx=None):
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

    if ctx is not None:
        ctx["magma_covar"] = f'{output_prefix}.gsa.out'
    print("\n🎉 MAGMA covariate analysis Completed.")
    return(f'{output_prefix}.gsa.out')