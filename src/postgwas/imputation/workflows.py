import argparse
import sys
import os
from pathlib import Path
from rich_argparse import RichHelpFormatter





# =========================================================
# DISPATCHERS
# =========================================================
def run_sumstat_imputation_direct(args):
    """
    DIRECT mode:
      - Run ONLY the chosen imputation tool.
      - No mandatory post-processing.
    """
    tool = args.imputation_tool
    if tool == "pred_ld":
        # Lazy import: only bring in PRED-LD when needed
        from postgwas.imputation.pred_ld.pred_ld_runner import (
            run_pred_ld_parallel,process_pred_ld_results_all_parallel
        )
        print("🚀 Running PRED-LD in DIRECT mode…")
        run_pred_ld_parallel(
            predld_input_dir=args.predld_input_dir,
            pred_ld_ref=args.ref_ld,
            output_folder=args.outdir,
            output_prefix=args.sample_id,
            r2threshold=args.r2threshold,
            maf=args.maf,
            population=args.population,
            ref=args.ref,
            threads=args.nthreads
        )
        print("✅ PRED-LD DIRECT run finished.")
    else:
        raise ValueError(f"Imputation tool not yet implemented for DIRECT mode: {tool}")
    process_pred_ld_results_all_parallel(
        folder_path=args.outdir,
        output_path=args.outdir,
        gwas2vcf_resource_folder=args.gwas2vcf_resource,
        output_prefix=args.sample_id,
        corr_method= "pearson",
        threads=args.nthreads,
    )

def run_sumstat_imputation_pipeline(args):
    """
    PIPELINE mode:
      - Run imputation
      - PLUS post-processing (merge + per-chromosome correlations).
    """
    tool = args.imputation_tool
    if tool == "pred_ld":
        from postgwas.imputation.pred_ld.pred_ld_runner import (
            run_pred_ld_parallel,
        )
        # from postgwas.imputation.pred_ld.pred_ld_result_processor import (
        #     process_pred_ld_results_all_parallel,
        # )
        print("🚀 Starting PRED-LD PIPELINE (imputation + post-processing)…")
        # 1) Run imputation
        run_pred_ld_parallel(
            predld_input_dir=args.predld_input_dir,
            pred_ld_ref=args.ref_ld,
            output_folder=args.outdir,
            sample_id=args.sample_id,
            r2threshold=args.r2threshold,
            maf=args.maf,
            population=args.population,
            ref=args.ref,
            threads=args.threads,
            verbose=args.verbose,
        )
        # 2) Post-processing: correlations + merging
        # print("\n📊 Running PRED-LD post-processing (correlations + merging)…")
        # combined_df, corr_df = process_pred_ld_results_all_parallel(
        #     folder_path=args.outdir,
        #     output_prefix=args.sample_id,
        #     corr_method=args.corr_method,
        #     threads=args.threads,
        # )
        # # 3) Optional: GWAS2VCF spec
        # if args.gwas2vcf_resource:
        #     print("\n🧩 Preparing GWAS2VCF specification file…")
        #     spec_path = os.path.join(
        #         args.outdir,
        #         f"{args.sample_id}_gwas2vcf_spec.tsv",
        #     )
        #     combined_df.to_csv(spec_path, sep="\t", index=False)
        #     print(f"💾 Saved GWAS2VCF spec file → {spec_path}")
        # print("\n✅ PRED-LD PIPELINE completed successfully!")
    else:
        raise ValueError(f"Imputation tool not yet implemented for PIPELINE mode: {tool}")
