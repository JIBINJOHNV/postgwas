#!/usr/bin/env python3
import argparse
from textwrap import dedent
from rich_argparse import RichHelpFormatter


class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    """Suppress default argparse error spam and show rich help instead."""
    def error(self, message):
        raise ValueError(message)


def build_finemap_parser():
    parser = HelpOnErrorArgumentParser(
        prog="finemap",
        formatter_class=RichHelpFormatter,
        description=dedent("""
        [bold cyan]FINEMAP v1.4.2[/bold cyan]
        [dim](c) 2015–2022 University of Helsinki[/dim]

        [bold]Help[/bold]
          • ./finemap --help
          • www.finemap.me
          • www.christianbenner.com

        [bold]Contact[/bold]
          • finemap@christianbenner.com
          • matti.pirinen@helsinki.fi
        """),
    )

    # ===============================================================
    # Subprograms (FINEMAP modes)
    # ===============================================================
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--sss",
        action="store_true",
        help="Fine-mapping with [bold]shotgun stochastic search[/bold]"
    )
    mode.add_argument(
        "--cond",
        action="store_true",
        help="Fine-mapping with [bold]stepwise conditional search[/bold]"
    )
    mode.add_argument(
        "--config",
        action="store_true",
        help="Evaluate a single causal configuration"
    )

    # ===============================================================
    # Input / execution
    # ===============================================================
    parser.add_argument(
        "--in-files",
        metavar="MASTER_FILE",
        required=True,
        help=(
            "Semicolon-separated master file with columns: "
            "[italic]'z', 'ld', 'snp', 'config', 'n_samples'[/italic] "
            "and optionally [italic]'k', 'log'[/italic]"
        )
    )

    parser.add_argument(
        "--dataset",
        default="all",
        help="Dataset selection (e.g. [italic]1,2[/italic] or [italic]1|2[/italic])"
    )

    parser.add_argument(
        "--n-threads",
        type=int,
        default=1,
        help="Number of parallel threads"
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Write output to log files defined in master file"
    )

    # ===============================================================
    # Shotgun stochastic search (SSS)
    # ===============================================================
    parser.add_argument(
        "--n-iter",
        type=int,
        default=100000,
        help="Maximum number of SSS iterations"
    )

    parser.add_argument(
        "--n-conv-sss",
        type=int,
        default=1000,
        help="Iterations used to check SSS convergence"
    )

    parser.add_argument(
        "--prob-conv-sss-tol",
        type=float,
        default=0.001,
        help="SSS convergence tolerance on added probability mass"
    )

    parser.add_argument(
        "--n-configs-top",
        type=int,
        default=50000,
        help="Number of top causal configurations to save"
    )

    # ===============================================================
    # Causal priors
    # ===============================================================
    parser.add_argument(
        "--n-causal-snps",
        type=int,
        default=5,
        help="Maximum number of causal SNPs"
    )

    parser.add_argument(
        "--prior-k",
        action="store_true",
        help="Use prior probabilities for number of causal SNPs from K files"
    )

    parser.add_argument(
        "--prior-k0",
        type=float,
        default=0.0,
        help="Prior probability that there are zero causal SNPs"
    )

    parser.add_argument(
        "--prior-snps",
        action="store_true",
        help="Use per-SNP causal prior probabilities from Z file"
    )

    parser.add_argument(
        "--prior-std",
        default="0.05",
        help="Comma-separated prior SDs of effect sizes"
    )

    # ===============================================================
    # LD / correlation handling
    # ===============================================================
    parser.add_argument(
        "--corr-config",
        type=float,
        default=0.95,
        help="Zero posterior if any SNP pair exceeds this absolute correlation"
    )

    parser.add_argument(
        "--collinear-tol",
        type=float,
        default=0.95,
        help="Collinearity tolerance for conditional search"
    )

    parser.add_argument(
        "--force-n-samples",
        action="store_true",
        help="Allow LD computed using different sample size than GWAS"
    )

    # ===============================================================
    # Filtering / credibility
    # ===============================================================
    parser.add_argument(
        "--pvalue-snps",
        type=float,
        default=1.0,
        help="P-value threshold for SNP inclusion"
    )

    parser.add_argument(
        "--cond-pvalue",
        type=float,
        default=5e-8,
        help="Genome-wide significance threshold (conditional search)"
    )

    parser.add_argument(
        "--prob-cred-set",
        type=float,
        default=0.95,
        help="Probability defining credible sets"
    )

    # ===============================================================
    # Configuration mode
    # ===============================================================
    parser.add_argument(
        "--rsids",
        help="Comma-separated list of rsIDs ([bold]required with --config[/bold])"
    )

    # ===============================================================
    # Misc
    # ===============================================================
    parser.add_argument(
        "--flip-beta",
        action="store_true",
        help="Flip effect sizes using 'flip' column in Z file"
    )

    parser.add_argument(
        "--std-effects",
        action="store_true",
        help="Print posterior mean and SD of standardized effect sizes"
    )

    return parser


# ======================================================================
# Entry point
# ======================================================================
if __name__ == "__main__":
    try:
        parser = build_finemap_parser()
        args = parser.parse_args()

        if args.config and not args.rsids:
            parser.error("--rsids is required when using --config")

        print("\n[✓] FINEMAP arguments parsed successfully\n")
        print(args)

    except ValueError:
        print("\n[bold red]❌ Invalid FINEMAP arguments[/bold red]\n")
        build_finemap_parser().print_help()
        raise SystemExit(2)
