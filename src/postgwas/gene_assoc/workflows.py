from pathlib import Path

from postgwas.formatter.cli import get_formatter_parser
from postgwas.gene_assoc.magma_main import magma_analysis_pipeline





# =====================================================================
#   DIRECT MODE ENGINE
# =====================================================================
def run_magma_direct(args):
    print("\n=== Running MAGMA Gene & Geneset Analysis (Direct Mode) ===")

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    magma_analysis_pipeline(
        output_dir=args.outdir,
        sample_id=args.sample_id,
        ld_ref=args.ld_ref,
        gene_loc_file=args.gene_loc_file,
        snp_loc_file=args.snp_loc_file,
        pval_file=args.pval_file,
        geneset_file=args.geneset_file,
        log_file=f"{args.outdir}/{args.sample_id}_magma_gene_assoc.log",
        threads=args.nthreads,
        num_batches=args.num_batches,
        window_upstream=args.window_upstream,
        window_downstream=args.window_downstream,
        gene_model=args.gene_model,
        n_sample_col=args.n_sample_col,
        seed=args.seed,
        magma=args.magma,
    )

    print("\n🎉 MAGMA Direct Mode Completed.")
    print(f"Results saved in: {args.outdir}\n")




# # =====================================================================
# #   FULL PIPELINE MODE (Step 1 → Step 5)
# # =====================================================================
# def run_magma_pipeline(args):
#     print("\n🚀 Running Full PostGWAS Pipeline → Harmonisation → LD Blocks → QC → Formatter → MAGMA")

#     # Step 1 + 2
#     run_annot_ldblock_pipeline(args)

#     # Step 3
#     run_sumstat_qc_pipeline(args)

#     # Step 4
#     run_formatter_pipeline(args)

#     # Step 5: MAGMA
#     print("\n--- Running MAGMA Step (5/5) ---")
#     Path(args.outdir).mkdir(parents=True, exist_ok=True)

#     magma_analysis_pipeline(
#         output_dir=args.outdir,
#         sample_id=args.sample_id,
#         ld_ref=args.ld_ref,
#         gene_loc_file=args.gene_loc_file,
#         snp_loc_file=args.snp_loc_file,
#         pval_file=args.pval_file,
#         geneset_file=args.geneset_file,
#         log_file=args.log_file if hasattr(args, "log_file") else None,
#         threads=args.nthreads,
#         num_batches=args.num_batches,
#         window_upstream=args.window_upstream,
#         window_downstream=args.window_downstream,
#         gene_model=args.gene_model,
#         n_sample_col=args.n_sample_col,
#         seed=args.seed,
#         magma=args.magma,
#     )

#     print("\n🎉 Full MAGMA Pipeline Completed.")
#     print(f"All results saved in: {args.outdir}\n")

