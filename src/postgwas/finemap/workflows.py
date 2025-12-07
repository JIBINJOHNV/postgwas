import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path

# SuSiE backend (Python) ---------------------
from postgwas.finemap.main import (
    validate_locus_file,
    run_susie,
)






# ===========================================
#   DIRECT MODE ENGINE
# ===========================================
def run_susie_direct(args):
    print("\n=== Running SuSiE Fine-Mapping (Direct Mode) ===")
    # Validate locus file
    validate_locus_file(args.locus_file)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    output_prefix = Path(args.outdir)

    run_susie(
        locus_file=args.locus_file,
        sumstat_file=args.finemap_input_file,
        sample_id=args.sample_id,
        ld_ref=args.finemap_ld_ref,
        plink=args.plink,
        output_folder=str(output_prefix),
        lp_threshold=args.lp_threshold,
        L=args.L,
        workers=args.nthreads,
        min_ram_per_worker_gb=args.min_ram_per_worker_gb,
        timeout_ld_seconds=args.timeout_ld_seconds,
        timeout_susie_seconds=args.timeout_susie_seconds,
        skip_mhc=args.finemap_skip_mhc,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
    )

    print("\n🎉 SuSiE Direct Mode Completed.")
    print(f"Results written to: {output_prefix}\n")


# ===========================================
#   PIPELINE MODE ENGINE
# ===========================================
def _infer_locus_file(args):
    """
    Decide how locus files are obtained in pipeline mode.
    Options:
      1) User-provided
      2) Generated during formatter / annotation step
    """
    if getattr(args, "locus_file", None):
        return args.locus_file

    # Default: formatter creates locus files inside the outdir
    candidate = Path(args.outdir) / f"{args.sample_id}.loci.tsv"
    return str(candidate)


def _infer_sumstat_file(args):
    """
    SuSiE requires a sumstat file with:
        SNP, BETA, SE, P, CHR, BP, Alleles...
    The formatter pipeline produces this file.
    """
    candidate = Path(args.outdir) / f"{args.sample_id}.susie_sumstats.tsv.gz"
    return str(candidate)





def run_susie_pipeline(args):

    print("\n🚀 Running Full PostGWAS SuSiE Pipeline → Harmonisation → LD → QC → Formatter → SuSiE")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Step 1 --------------------------------
    print("\n[1/5] Harmonisation…")
    run_harmonisation_pipeline(args)

    # Step 2 --------------------------------
    print("\n[2/5] LD Block Annotation…")
    run_annot_ldblock_pipeline(args)

    # Step 3 --------------------------------
    print("\n[3/5] Summary Statistics QC…")
    run_sumstat_qc_pipeline(args)

    # Step 4 --------------------------------
    print("\n[4/5] Formatter → SuSiE-ready inputs…")
    run_formatter_pipeline(args)

    # Prepare SuSiE inputs -------------------
    locus_file = _infer_locus_file(args)
    sumstat_file = _infer_sumstat_file(args)

    print("\n[5/5] SuSiE Fine-Mapping…")
    validate_locus_file(locus_file)

    output_prefix = out / args.sample_id

    run_susie(
        locus_file=locus_file,
        sumstat_file=sumstat_file,
        sample_id=args.sample_id,
        ld_ref=args.finemap_ld_ref,
        plink=args.plink,
        output_folder=str(output_prefix),
        lp_threshold=args.lp_threshold,
        L=args.L,
        workers=args.nthreads,
        min_ram_per_worker_gb=args.min_ram_per_worker_gb,
        timeout_ld_seconds=args.timeout_ld_seconds,
        timeout_susie_seconds=args.timeout_susie_seconds,
        skip_mhc=args.skip_mhc,
        mhc_start=args.mhc_start,
        mhc_end=args.mhc_end,
    )

    print("\n🎉 SuSiE Pipeline Completed.")
    print(f"All results saved under: {output_prefix}\n")




def run_finemap_direct(args):
    pass

def run_finemap_pipeline(args):
    pass 