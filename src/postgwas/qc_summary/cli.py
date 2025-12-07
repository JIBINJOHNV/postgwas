#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------
from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_common_out_parser,
    get_bcftools_binary_parser,
)

from postgwas.harmonisation.cli import get_harmonisation_parser
from postgwas.annot_ldblock.cli import get_annot_ldblock_parser



# Import your QC summary function
from postgwas.qc_summary.workflows import run_qc_summary_direct,run_qc_summary_pipeline



# ==================================================================
#  MAIN CLI
# ==================================================================
def main():
    top_parser = argparse.ArgumentParser(
        prog="qc_summary",
        description=(
            "QC Summary Module\n\n"
            "This module performs summary-level QC on VCF files using bcftools. "
            "It extracts essential variant-quality metrics using:\n"
            "  • bcftools view (AF checks)\n"
            "  • bcftools stats\n\n"
            "It can run in:\n"
            "  • DIRECT mode — QC on an existing harmonised/filtered VCF\n"
            "  • PIPELINE mode — Run Steps 1–3 (Harmonisation → LDBlock → SumstatFilter) "
            "    and then run QC Summary.\n"
        ),
        formatter_class=RichHelpFormatter
    )

    subparsers = top_parser.add_subparsers(dest="mode", help="direct or pipeline")

    # ---------------------------------------------------------
    # DIRECT MODE
    # ---------------------------------------------------------
    direct_parser = subparsers.add_parser(
        "direct",
        help="Run QC summary on an existing merged/filtered VCF.",
        parents=[
            get_defaultresourse_parser(),
            get_inputvcf_parser(),
            get_common_out_parser(),
            get_bcftools_binary_parser()
        ],
        formatter_class=RichHelpFormatter
    )

    direct_parser.add_argument(
        "--external-af-name",
        metavar=" ",
        default="EUR",
        help=(
            "Name of the INFO/<AF> tag used as external reference AF.\n"
            "Examples: EUR, AFR, EAS.\n"
            "Default: EUR"
        )
    )

    direct_parser.add_argument(
        "--allelefreq-diff-cutoff",
        type=float,
        metavar=" ",
        default=0.2,
        help=(
            "Maximum allowed absolute difference between FORMAT/AF and INFO/<external_af_name>.\n"
            "Default: 0.2"
        )
    )

    direct_parser.set_defaults(func=run_qc_summary_direct)

    # ---------------------------------------------------------
    # PIPELINE MODE
    # ---------------------------------------------------------
    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Run Harmonisation → LD Block Annot → Sumstat Filter → QC Summary.",
        parents=[
            get_defaultresourse_parser(),
            get_harmonisation_parser(),
            get_annot_ldblock_parser(),
            get_common_out_parser(),
            get_bcftools_binary_parser()
        ],
        formatter_class=RichHelpFormatter
    )
    pipeline_parser.set_defaults(func=run_qc_summary_pipeline)

    # ---------------------------------------------------------
    # EXECUTE
    # ---------------------------------------------------------
    args = top_parser.parse_args()

    # Avoid running with no args → show help
    if len(sys.argv) == 1:
        top_parser.print_help()
        sys.exit(0)

    if args.mode == "direct" and len(sys.argv) == 2:
        top_parser.parse_args(["direct", "--help"])
        sys.exit(0)

    if args.mode == "pipeline" and len(sys.argv) == 2:
        top_parser.parse_args(["pipeline", "--help"])
        sys.exit(0)

    # Execute chosen mode
    if hasattr(args, "func"):
        args.func(args)
    else:
        top_parser.print_help()


if __name__ == "__main__":
    main()
