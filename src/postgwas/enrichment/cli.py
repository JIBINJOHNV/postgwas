#!/usr/bin/env python3
import argparse
from postgwas.utils.main import validate_path
from rich_argparse import RichHelpFormatter
from postgwas.clis.geneset_cli import get_geneset_common_parser,get_common_out_parser, get_defaultresourse_parser


# =========================================================
# SHARED LDSC INPUT ARGUMENTS
# =========================================================
def get_geneset_parser(add_help: bool = False) -> argparse.ArgumentParser:
    """
    Defines input arguments for gene set over representation analysis.
    """
    parser = argparse.ArgumentParser(
        add_help=add_help,
        formatter_class=RichHelpFormatter,
    )

    # Attach shared argument groups
    get_geneset_common_parser()
    get_common_out_parser()
    get_defaultresourse_parser()
    return parser


