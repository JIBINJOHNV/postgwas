
from postgwas.sumstat_filter.sumstat_filter import filter_gwas_vcf_bcftools
from pathlib import Path


def run_sumstat_filter_direct(args):
    # Your logic — YOU said do not enforce checks here
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    filter_gwas_vcf_bcftools(
        vcf_path=args.vcf,
        output_folder=str(outdir),
        output_prefix=args.sample_id,
        pval_cutoff=args.pval_cutoff,
        maf_cutoff=args.maf_cutoff,
        allelefreq_diff_cutoff=args.allelefreq_diff_cutoff,
        info_cutoff=args.info_cutoff,
        external_af_name=args.external_af_name,
        include_indels=args.include_indels,
        include_palindromic=args.include_palindromic,
        remove_mhc=args.remove_mhc,
        mhc_chrom=args.mhc_chrom,
        mhc_start=args.mhc_start,
        mhc_end=args.mhc_end,
        threads=args.nthreads,
        max_mem=args.max_mem
    )

    print("✅ Step 3 QC Filtering Completed.")


# ---------------------------------------------------------
# Step-1 → Step-2 → Step-3 pipeline
# ---------------------------------------------------------
def run_sumstat_filter_pipeline(args):
    print("🚀 Running Harmonisation → Annot LDBlock → Sumstat QC Pipeline")
#     run_annot_ldblock_pipeline(args)  # includes step-1 + step-2
#     run_sumstat_qc(args)
#     print("🎉 Full Pipeline Completed.")
