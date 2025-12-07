#!/usr/bin/env python3

from pathlib import Path
import os
import subprocess
from postgwas.formatter.to_magma import vcf_to_magma
from postgwas.gene_tests.magma_main import magma_analysis_pipeline
from postgwas.magmacovar.main import run_magma_covariates
from postgwas.finemap.cli import main as finemap_main
from postgwas.pops.pops import pops_main, get_pops_args
from postgwas.finemapping.main import validate_locus_file, run_susie


# vcf_to_magma(sumstat_vcf=, 
#              output_folder=
#              sample_name=)

# magma_analysis_pipeline(
#     output_dir=,
#     sample_id=,
#     ld_ref=,
#     gene_loc_file=,
#     snp_loc_file=,
#     pval_file=,
#     geneset_file=,
#     log_file=,
#     threads= 6,
#     num_batches= 6,
#     window_upstream= 35,
#     window_downstream= 10,
#     gene_model= "snp-wise=mean",
#     n_sample_col= "N_COL",
#     seed= 10,
#     magma= "magma"
# )

# run_magma_covariates(
#     magma_bin = ,
#     gene_results_file =,
#     covariates_file =,
#     output_prefix =,
#     model= None,
#     direction= "two-sided",
#     log_file= None,
# )

# run_susie(
#     locus_file,
#     sumstat_file,
#     sample_id,
#     ld_ref,
#     plink,
#     output_folder,
#     lp_threshold="7.3",
#     L="10",
#     workers="auto",
#     min_ram_per_worker_gb="4",
#     timeout_ld_seconds="180",
#     timeout_susie_seconds="180",
#     skip_mhc=False,
#     mhc_start="25000000",
#     mhc_end="35000000",
#     verbose=False,
# )




'''

    # Construct PoPS argument list
    pops_argv = [
        "--gene_annot_path", args.pops_gene_loc_file,
        "--feature_mat_prefix", args.feature_mat_prefix,
        "--num_feature_chunks", str(args.num_feature_chunks),
        "--magma_prefix", magma_output["merged_prefix"],
        "--out_prefix", f'{pops_folder}/{args.sample_name}_pops',
    ]

    if args.control_features_path:
        pops_argv += ["--control_features_path", args.control_features_path]

    pops_argv.append("--use_magma_covariates" if args.use_magma_covariates else "--ignore_magma_covariates")
    pops_argv.append("--use_magma_error_cov" if args.use_magma_error_cov else "--ignore_magma_error_cov")

    if args.y_path:
        pops_argv += ["--y_path", args.y_path]
    if args.y_covariates_path:
        pops_argv += ["--y_covariates_path", args.y_covariates_path]
    if args.y_error_cov_path:
        pops_argv += ["--y_error_cov_path", args.y_error_cov_path]

    if args.project_out_covariates_chromosomes:
        pops_argv += ["--project_out_covariates_chromosomes"] + args.project_out_covariates_chromosomes

    pops_argv.append("--project_out_covariates_remove_hla" if args.project_out_covariates_remove_hla else "--project_out_covariates_keep_hla")

    if args.subset_features_path:
        pops_argv += ["--subset_features_path", args.subset_features_path]

    if args.feature_selection_chromosomes:
        pops_argv += ["--feature_selection_chromosomes"] + args.feature_selection_chromosomes

    pops_argv += ["--feature_selection_p_cutoff", str(args.feature_selection_p_cutoff)]

    if args.feature_selection_max_num is not None:
        pops_argv += ["--feature_selection_max_num", str(args.feature_selection_max_num)]

    if args.feature_selection_fss_num_features is not None:
        pops_argv += ["--feature_selection_fss_num_features", str(args.feature_selection_fss_num_features)]

    pops_argv.append("--feature_selection_remove_hla" if args.feature_selection_remove_hla else "--feature_selection_keep_hla")

    if args.training_chromosomes:
        pops_argv += ["--training_chromosomes"] + args.training_chromosomes

    pops_argv.append("--training_remove_hla" if args.training_remove_hla else "--training_keep_hla")

    pops_argv += ["--method", args.method]
    pops_argv.append("--save_matrix_files" if args.save_matrix_files else "--no_save_matrix_files")
    pops_argv += ["--random_seed", str(args.random_seed)]
    pops_argv.append("--verbose" if args.verbose else "--no_verbose")

    # Pass args to PoPS
    pops_args = get_pops_args(pops_argv)
    return pops_main(vars(pops_args))

'''





def run_flames(
    annotation_dir,
    pops,
    magma_z,
    magma_tissue,
    indexfile,
    outdir,
    filename="FLAMES_scores",
    modelpath=None,
    snp_col="cred1",
    prob_col="prob1",
    distance=750000,
    weight=0.725,
):

    # Auto-detect FLAMES.py path
    flames_py = Path(__file__).parent / "FLAMES.py"

    # -----------------------------
    # 🔥 STEP 1 — FLAMES annotate
    # -----------------------------
    cmd_annot = f"""
    python "{flames_py}" annotate \
        --annotation_dir "{annotation_dir}" \
        --pops "{pops}" \
        --magma_z "{magma_z}" \
        --magma_tissue "{magma_tissue}" \
        --indexfile "{indexfile}" \
        --prob_col "{prob_col}" \
        --SNP_col "{snp_col}"
    """

    print("🔥 Running FLAMES annotation…")
    subprocess.run(cmd_annot, shell=True, check=True)

    # -----------------------------
    # 🔥 STEP 2 — FLAMES scoring
    # -----------------------------
    cmd_score = f"""
    python "{flames_py}" FLAMES \
        --indexfile "{indexfile}" \
        --outdir "{outdir}" \
        --filename "{filename}" \
        --distance {distance} \
        --weight {weight} \
        --modelpath "{modelpath if modelpath else Path(__file__).parent / 'model'}"
    """

    print("🔥 Running FLAMES scoring…")
    subprocess.run(cmd_score, shell=True, check=True)

    print("🎉 FLAMES annotation + scoring completed.")


if __name__ == "__main__":





# def run_pipeline(args):
#     outdir = Path(args.out)
#     outdir.mkdir(parents=True, exist_ok=True)

#     # ============================================================
#     # MODE 1: DIRECT → Only FLAMES Annotate + FLAMES
#     # ============================================================
#     if args.mode == "direct":
#         print("🔥 DIRECT mode: Running only FLAMES annotate → FLAMES")

#         functional_annotation(
#             annotation_dir=args.flames_annot,
#             magma_genes=args.magma_genes,
#             magma_ts=args.magma_ts,
#             pops=args.pops,
#             finemap=args.finemap_cred,
#             outdir=outdir / "flames_annotation"
#         )

#         FLAMES(
#             annotation_dir=outdir / "flames_annotation",
#             outdir=outdir / "flames_final"
#         )

#         print("\n🎉 DIRECT mode completed.")
#         return

#     # ============================================================
#     # MODE 2: PIPELINE → Full process starting from VCF
#     # ============================================================
#     if args.mode == "pipeline":
#         print("🚀 PIPELINE mode: Full workflow starting from VCF")

#         # ------------------ Step 1: VCF to MAGMA/FineMap Inputs ------------------
#         print("\n📌 Step 1: Processing VCF…")
#         inputs = vcf_to_magma(
#             vcf=args.vcf,
#             outdir=outdir / "inputs",
#             gene_loc=args.gene_loc,
#             resource_folder=args.resource_folder
#         )
#         snp_loc = inputs["snp_loc_file"]
#         pval_file = inputs["pval_file"]

#         # ------------------ Step 2: MAGMA ------------------
#         print("\n📌 Step 2: Running MAGMA gene analysis…")
#         magma_genes = run_magma(
#             outdir / "magma",
#             inputs=snp_loc,
#             pval=pval_file,
#             ld_ref=args.ld_panel,
#             gene_loc=args.gene_loc
#         )

#         # ------------------ Step 3: MAGMA Tissue ------------------
#         print("\n📌 Step 3: Running MAGMA tissue…")
#         magma_ts = run_magma_tissue(
#             outdir / "magma_tissue",
#             magma_genes,
#             resource_folder=args.resource_folder
#         )

#         # ------------------ Step 4: FineMap ------------------
#         print("\n📌 Step 4: Running fine-mapping…")
#         finemap_cred = run_finemap(
#             outdir / "finemap",
#             vcf=args.vcf,
#             ld_panel=args.ld_panel
#         )

#         # ------------------ Step 5: PoPS ------------------
#         print("\n📌 Step 5: Running PoPS…")
#         pops_preds = run_pops(
#             outdir / "pops",
#             pval_file=pval_file,
#             resource_folder=args.resource_folder
#         )

#         # ------------------ Step 6: FLAMES Annotate ------------------
#         print("\n📌 Step 6: Running FLAMES annotate…")
#         flames_annotate(
#             annotation_dir=args.flames_annot,
#             magma_genes=magma_genes,
#             magma_ts=magma_ts,
#             pops=pops_preds,
#             finemap=finemap_cred,
#             outdir=outdir / "flames_annotation"
#         )

#         # ------------------ Step 7: FLAMES ------------------
#         print("\n📌 Step 7: Running FLAMES…")
#         flames_run(
#             annotation_dir=outdir / "flames_annotation",
#             outdir=outdir / "flames_final"
#         )

#         print("\n🎉 PIPELINE mode completed.")
