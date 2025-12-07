from postgwas.qc_summary.main import run_qc_summary



def run_qc_summary_direct(args):
    run_qc_summary(
        vcf_path=args.vcf,
        qc_outdir=args.outdir,
        sample_id=args.sample_id,
        external_af_name=args.external_af_name,
        allelefreq_diff_cutoff=args.allelefreq_diff_cutoff,
        n_threads=args.nthreads,
        bcftools_bin=args.bcftools,
    )





# ==================================================================
#  PIPELINE MODE FUNCTION
# ==================================================================
def run_qc_summary_pipeline(args):
    """
    Pipeline QC Summary:
    1. Harmonisation
    2. LD-Block Annotation
    3. Sumstat Filtering
    4. QC Summary

    This mirrors the structure of other PostGWAS pipeline CLIs.
    """

    # ---------------------------------
    # Step 1 — Harmonisation
    # ---------------------------------
    harm_df = run_harmonisation(args)

    # ---------------------------------
    # Step 2 — Annotate LD Blocks
    # ---------------------------------
    run_annot_ldblock_pipeline(args)

    # ---------------------------------
    # Step 3 — Sumstat Filter
    # ---------------------------------
    filtered_vcf, _, _ = run_sumstat_filter_pipeline(args)

    # ---------------------------------
    # Step 4 — QC Summary
    # ---------------------------------
    run_qc_summary(
        vcf_path=filtered_vcf,
        qc_outdir=args.outdir,
        sample_id=args.sample_id,
        external_af_name=args.external_af_name,
        allelefreq_diff_cutoff=args.allelefreq_diff_cutoff,
        n_threads=args.nthreads,
        bcftools_bin=args.bcftools,
    )

