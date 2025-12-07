
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---------------------------------------------------------
# Shared CLI components (must accept add_help=False)
# ---------------------------------------------------------
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_genomeversion_parser,
    get_common_out_parser,
    get_common_magma_assoc_parser,
    get_common_pops_parser,     # This contains all PoPS arguments
)

# Pipeline step parsers
from postgwas.harmonisation.cli import (
    get_harmonisation_parser,
)
from postgwas.annot_ldblock.cli import (
    get_annot_ldblock_parser,
)




import sys
import argparse
from pathlib import Path

# Internal imports
from postgwas.formatter.to_magma import vcf_to_magma
from postgwas.gene_assoc.magma_main import run_magma_analysis
from postgwas.pops.pops import pops_main, get_pops_args



def run_pops_direct(args: argparse.Namespace) -> None:
    """
    Run PoPS directly using provided arguments.
    This bypasses the full PostGWAS pipeline and runs PoPS only.
    """
    print("\n📌 PoPS — Direct Mode")
    # --------------------------------------------
    # Create output folder
    # --------------------------------------------
    pops_folder = Path(args.outdir)
    pops_folder.mkdir(parents=True, exist_ok=True)

    pops_out_prefix = pops_folder / f"{args.sample_id}_pops"

    # --------------------------------------------
    # Build argument list for PoPS
    # --------------------------------------------
    pops_argv = [
        "--gene_annot_path", args.pops_gene_loc_file,
        "--feature_mat_prefix", args.feature_mat_prefix,
        "--num_feature_chunks", str(args.num_feature_chunks),
        "--magma_prefix", args.magma_assoc_prefix,
        "--out_prefix", str(pops_out_prefix),
    ]

    # Optional: Control features
    if getattr(args, "control_features_path", None):
        pops_argv += ["--control_features_path", args.control_features_path]

    # Optional flags: use or ignore MAGMA covariates
    pops_argv.append(
        "--use_magma_covariates" if args.use_magma_covariates else "--ignore_magma_covariates"
    )
    pops_argv.append(
        "--use_magma_error_cov" if args.use_magma_error_cov else "--ignore_magma_error_cov"
    )

    # Optional Y data
    if getattr(args, "y_path", None):
        pops_argv += ["--y_path", args.y_path]
    if getattr(args, "y_covariates_path", None):
        pops_argv += ["--y_covariates_path", args.y_covariates_path]
    if getattr(args, "y_error_cov_path", None):
        pops_argv += ["--y_error_cov_path", args.y_error_cov_path]

    # Project-out chromosomes
    if getattr(args, "project_out_covariates_chromosomes", None):
        pops_argv += [
            "--project_out_covariates_chromosomes",
            *args.project_out_covariates_chromosomes
        ]

    pops_argv.append(
        "--project_out_covariates_remove_hla"
        if args.project_out_covariates_remove_hla
        else "--project_out_covariates_keep_hla"
    )

    # Subset features
    if getattr(args, "subset_features_path", None):
        pops_argv += ["--subset_features_path", args.subset_features_path]

    # Feature selection chromosomes
    if getattr(args, "feature_selection_chromosomes", None):
        pops_argv += [
            "--feature_selection_chromosomes",
            *args.feature_selection_chromosomes,
        ]

    pops_argv += [
        "--feature_selection_p_cutoff",
        str(args.feature_selection_p_cutoff)
    ]

    if args.feature_selection_max_num is not None:
        pops_argv += [
            "--feature_selection_max_num",
            str(args.feature_selection_max_num),
        ]

    if args.feature_selection_fss_num_features is not None:
        pops_argv += [
            "--feature_selection_fss_num_features",
            str(args.feature_selection_fss_num_features),
        ]

    pops_argv.append(
        "--feature_selection_remove_hla"
        if args.feature_selection_remove_hla
        else "--feature_selection_keep_hla"
    )

    if getattr(args, "training_chromosomes", None):
        pops_argv += [
            "--training_chromosomes",
            *args.training_chromosomes,
        ]

    pops_argv.append(
        "--training_remove_hla"
        if args.training_remove_hla
        else "--training_keep_hla"
    )

    pops_argv += [
        "--method", args.method,
        "--random_seed", str(args.seed),
    ]

    pops_argv.append(
        "--save_matrix_files"
        if args.save_matrix_files
        else "--no_save_matrix_files"
    )

    pops_argv.append(
        "--verbose"
        if args.pops_verbose
        else "--no_verbose"
    )

    # --------------------------------------------
    # Run PoPS
    # --------------------------------------------
    print(f"\n🚀 Running PoPS with {len(pops_argv)} arguments…")
    pops_args = get_pops_args(pops_argv)
    pops_main(vars(pops_args))

    print("\n🎉 PoPS Direct Mode Completed Successfully!")
    print(f"Output prefix: {pops_out_prefix}\n")



def run_pops_pipeline(args):
    """
    Executes the VCF -> MAGMA -> PoPS workflow.
    args: Namespace object containing all necessary parameters.
    """
    print("\n====================")
    print("      PIPELINE MODE")
    print("====================")

    # ---------- Step 1: VCF to MAGMA ----------
    print("\n📌 Step 1 — Preparing MAGMA inputs from VCF...")
    magma_inputs = vcf_to_magma(
        sumstat_vcf=args.sumstat_vcf,
        output_folder=f"{args.outdir}/magma_inputs/",
        sample_id=args.sample_id,
    )

    snp_loc_file = magma_inputs["snp_loc_file"]
    pval_file = magma_inputs["pval_file"]

    # Create MAGMA analysis folder
    magma_folder = f"{args.outdir}/magma_analysis/"
    Path(magma_folder).mkdir(parents=True, exist_ok=True)

    # ---------- Step 2: Run MAGMA ----------
    print("\n📌 Step 2 — Running MAGMA...")
    log_file = str(Path(magma_folder) / f"{args.sample_id}.magma.log")

    magma_output = run_magma_analysis(
        magma_analysis_folder=magma_folder,
        sample_id=args.sample_id,
        ld_ref=args.ld_ref,
        gene_loc_file=args.magma_gene_loc_file,
        snp_loc_file=snp_loc_file,
        pval_file=pval_file,
        geneset_file=args.geneset_file,
        log_file=log_file,
        window_upstream=args.window_upstream,
        window_downstream=args.window_downstream,
        gene_model=args.gene_model,
        n_sample_col=args.n_sample_col,
        num_batches=args.num_batches,
        num_cores=args.num_cores,
        seed=args.seed,
        magma=args.magma,
    )

    # ---------- Step 3: Run PoPS ----------
    print("\n📌 Step 3 — Running PoPS...")

    pops_folder = f"{args.outdir}/"
    Path(pops_folder).mkdir(parents=True, exist_ok=True)
    
    # Construct PoPS argument list. instead of this arg need to update
    pops_argv = [
        "--gene_annot_path", args.pops_gene_loc_file,
        "--feature_mat_prefix", args.feature_mat_prefix,
        "--num_feature_chunks", str(args.num_feature_chunks),
        "--magma_prefix", magma_output["merged_prefix"],
        "--out_prefix", f'{pops_folder}/{args.sample_id}_pops',
    ]
    
    run_pops_direct(args)
