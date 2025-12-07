#!/usr/bin/env python3
"""
Python wrapper to call the Rscript CLI for SuSiE parallel fine-mapping.
"""

import argparse
import subprocess
from pathlib import Path
import sys
from postgwas.finemapping.main import validate_locus_file, run_susie


def main(argv=None):

    parser = argparse.ArgumentParser(
        prog="susie-finemap",
        description="Python wrapper for SuSiE fine-mapping (calls the R CLI script).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # -------------------------------------------------------------
    # Required arguments
    # -------------------------------------------------------------
    parser.add_argument("--locus_file", required=True)
    parser.add_argument("--sumstat_file", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--ld_ref", required=True)
    parser.add_argument("--plink", required=True)
    parser.add_argument("--SUSIE_Analysis_folder", required=True)

    # -------------------------------------------------------------
    # Optional arguments
    # -------------------------------------------------------------
    parser.add_argument("--lp_threshold", default="7.3")
    parser.add_argument("--L", default="10")
    parser.add_argument("--workers", default="auto")
    parser.add_argument("--min_ram_per_worker_gb", default="4")
    parser.add_argument("--timeout_ld_seconds", default="180")
    parser.add_argument("--timeout_susie_seconds", default="180")
    parser.add_argument("--skip_mhc", action="store_true")
    parser.add_argument("--mhc_start", default="25000000")
    parser.add_argument("--mhc_end", default="35000000")
    parser.add_argument("--verbose", action="store_true")

    # -------------------------------------------------------------
    # Parse args (CLI or programmatically)
    # -------------------------------------------------------------
    args = parser.parse_args(argv)

    # Validate locus file
    validate_locus_file(args.locus_file)

    # -------------------------------------------------------------
    # Run SuSiE (Python function wrapper)
    # -------------------------------------------------------------
    run_susie(
        locus_file=args.locus_file,
        sumstat_file=args.sumstat_file,
        sample_id=args.sample_id,
        ld_ref=args.ld_ref,
        plink=args.plink,
        output_folder=args.SUSIE_Analysis_folder,
        lp_threshold=args.lp_threshold,
        L=args.L,
        workers=args.workers,
        min_ram_per_worker_gb=args.min_ram_per_worker_gb,
        timeout_ld_seconds=args.timeout_ld_seconds,
        timeout_susie_seconds=args.timeout_susie_seconds,
        skip_mhc=args.skip_mhc,
        mhc_start=args.mhc_start,
        mhc_end=args.mhc_end,
        verbose=args.verbose
    )

    print("\n✅ SuSiE fine-mapping completed successfully.\n")


if __name__ == "__main__":
    main()
