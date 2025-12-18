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
    output_folder = Path(args.outdir)

    susie_results=run_susie(
        locus_file=args.locus_file,
        sumstat_file=args.finemap_input_file,
        sample_id=args.sample_id,
        ld_ref=args.finemap_ld_ref,
        plink=args.plink,
        output_folder=str(output_folder),
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
    print(f"\n🎉 SuSiE Direct Mode Completed. {susie_results}")
    return susie_results
