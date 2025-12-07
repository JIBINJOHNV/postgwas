#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---------------------------------------------------------
# IMPORT SHARED PARSERS
# ---------------------------------------------------------
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_bcftools_binary_parser,

)

from postgwas.harmonisation.cli import (
    get_harmonisation_parser,
    run_harmonisation,
)

from postgwas.annot_ldblock.cli import get_annot_ldblock_parser

from postgwas.ld_clump.workflows import (
    run_ld_clump_direct,
    run_ld_clump_pipeline,
)



def get_ld_clump_parser(add_help=False):
    """
    Returns a parent-safe parser containing arguments
    that are shared between DIRECT and PIPELINE modes.

    COMMON arguments:
        --ld-mode
        --population
        --r2-cutoff
        --window-kb

    (VCF input, outdir, resource folder, bcftools, harmonisation,
     annot_ldblock etc. are added separately in main() via their own
     get_*_parser functions.)
    """
    parser = argparse.ArgumentParser(add_help=add_help)

    # -------------------------------
    # General LD clumping options
    # -------------------------------
    general = parser.add_argument_group("Common LD Clumping Arguments")

    general.add_argument(
        "--ld-mode",
        default="by_regions",
        choices=["by_regions", "standard"],
        metavar=" ",
        help=(
            "LD clumping strategy.\n"
            "Choices: [cyan]by_regions[/cyan], [cyan]standard[/cyan]\n"
            f"[bold green]Default:[/bold green] [cyan]by_regions[/cyan]"
        )
    )

    general.add_argument(
        "--population",
        default="EUR",
        choices=["EUR", "AFR", "EAS"],
        metavar=" ",
        help=(
            "Population LD blocks (by_regions mode only).\n"
            f"[bold green]Default:[/bold green] [cyan]EUR[/cyan]"
        )
    )

    # -------------------------------
    # Standard-mode arguments
    # -------------------------------
    standard = parser.add_argument_group(
        "Standard LD Clumping Arguments (PLINK-style)"
    )

    standard.add_argument(
        "--r2-cutoff",
        type=float,
        default=0.1,
        metavar=" ",
        help=(
            "Pairwise LD r² threshold for standard clumping.\n"
            f"[bold green]Default:[/bold green] [cyan]0.1[/cyan]"
        )
    )

    standard.add_argument(
        "--window-kb",
        type=int,
        default=250,
        metavar=" ",
        help=(
            "Sliding window size (kilobases) for standard clumping.\n"
            f"[bold green]Default:[/bold green] [cyan]250[/cyan]"
        )
    )

    return parser




def main():
    # ---------------------------------------------------------
    # Top-level parser (Combined Direct/Pipeline arguments)
    # ---------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="ld_clump",
        description=(
            "LD Clumping Module (Direct Run Mode)\n\n"
            "Performs LD pruning using precomputed LD blocks or standard r² + window."
        ),
        formatter_class=RichHelpFormatter,
        parents=[
            # These are the combined parsers required for the full functionality:
            get_defaultresourse_parser(),
            get_bcftools_binary_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_ld_clump_parser(add_help=False)
        ]
    )

    # ---------------------------------------------------------
    # EXECUTE
    # ---------------------------------------------------------
    # If no arguments provided → show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    # Parse all arguments
    args = parser.parse_args()

    # We set the function to be executed explicitly:
    run_ld_clump_direct(args)


if __name__ == "__main__":
    # Ensure necessary imports are present for the script to run
    import argparse
    import sys
    main()





# ld_independent_snp(
#     sumstat_vcf="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_GRCh37_merged.vcf.gz",
#     output_folder="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/",
#     sample_name="PGC3_SCZ_european",
#     population="EUR"
# )