#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v10, Option B Error Handling)

Changes in this version
-----------------------
• argparse error messages fully suppressed
• custom message shown instead of argparse default:
      ❗ Required arguments are missing for the selected modules.
         Displaying full pipeline help:
• On missing args, pipeline plan + full help is shown dynamically
• Optional flags (--apply-filter, --apply-imputation, --apply-manhattan)
  behave as modules and expand dependencies automatically
"""

import argparse
import sys
from rich_argparse import RichHelpFormatter

# =====================================================================
# CUSTOM PARSER (SUPPRESS DEFAULT ARGPARSE ERROR OUTPUT)
# =====================================================================

class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        # We suppress default argparse error message,
        # print only our custom help block.
        raise ValueError("ARGPARSE_ERROR")  # handled later


# =====================================================================
# WORKFLOW IMPORTS
# =====================================================================

from postgwas.annot_ldblock.workflows import run_annot_ldblock
from postgwas.formatter.workflows import run_formatter_direct
from postgwas.sumstat_filter.workflows import run_sumstat_filter_direct
from postgwas.imputation.workflows import run_sumstat_imputation_direct
from postgwas.finemap.workflows import run_susie_direct
from postgwas.ld_clump.workflows import run_ld_clump_direct
from postgwas.gene_assoc.workflows import run_magma_direct
from postgwas.magmacovar.workflows import run_magma_covar_direct
from postgwas.pops.workflows import run_pops_direct
from postgwas.flames.workflows import run_flames_direct
from postgwas.h2_rg.workflows import run_ldsc_direct
from postgwas.manhattan.workflows import run_assoc_plot_direct
from postgwas.qc_summary.workflows import run_qc_summary_direct

# =====================================================================
# PARSER BUILDERS
# =====================================================================

from postgwas.clis.common_cli import (
    get_defaultresourse_parser,
    get_inputvcf_parser,
    get_genomeversion_parser,
    get_common_out_parser,
    get_common_sumstat_filter_parser,
    get_formatter_parser,
    get_finemap_method_parser,
    get_common_susie_arguments,
    get_plink_binary_parser,
    get_common_magma_assoc_parser,
    get_magma_binary_parser,
    get_common_magma_covar_parser,
    get_common_pops_parser,
    get_flames_common_parser,
    get_common_imputation_parser,
    get_ldsc_common_parser,
    get_annot_ldblock_parser,
    get_assoc_plot_parser
)

def get_qc_parser():
    parser = argparse.ArgumentParser(add_help=False)
    g = parser.add_argument_group("QC Summary Arguments")
    g.add_argument("--qc_output_prefix")
    return parser


# =====================================================================
# MODULE DESCRIPTIONS (for help preview)
# =====================================================================

MODULE_DESCRIPTIONS = {
    "annot_ldblock": [
        "Annotates variants with LD blocks using Berisa & Pickrell (2016).",
        "Inputs: harmonised VCF, genome build, LD block BED resources."
    ],
    "formatter": [
        "Converts harmonised VCF to tool-specific formats.",
        "Outputs: MAGMA, FINEMAP, LDSC, PoPS, FLAMES-ready tables."
    ],
    "sumstat_filter": [
        "Filters summary statistics: INFO, MAF, MAC, allele matching, EAF QC.",
        "Recommended before imputation and fine-mapping."
    ],
    "imputation": [
        "Performs summary statistics imputation using PRED-LD.",
        "Requires: LD reference panel matched to genome build."
    ],
    "finemap": [
        "Runs SuSiE fine-mapping to compute credible sets & PIPs.",
        "Requires: locus file, LD matrix, formatted summary statistics."
    ],
    "ld_clump": [
        "Runs PLINK LD-clumping to extract independent GWAS lead variants.",
    ],
    "magma": [
        "Runs MAGMA gene-level association analysis.",
    ],
    "magmacovar": [
        "Runs MAGMA gene-property (covariate) analysis.",
    ],
    "pops": [
        "Computes PoPS gene prioritisation scores."
    ],
    "flames": [
        "Runs FLAMES integrative effector-gene scoring."
    ],
    "heritability": [
        "Runs LDSC heritability & genetic correlation."
    ],
    "manhattan": [
        "Generates Manhattan & QQ plots."
    ],
    "qc": [
        "Generates harmonisation/QC summaries."
    ],
}

# =====================================================================
# DAG DEPENDENCIES
# =====================================================================

MODULE_DEPENDENCIES = {
    "annot_ldblock": [],
    "formatter": ["annot_ldblock"],
    "sumstat_filter": ["annot_ldblock"],
    "imputation": ["sumstat_filter"],
    "finemap": ["formatter", "annot_ldblock"],
    "ld_clump": ["formatter"],
    "magma": ["formatter"],
    "magmacovar": ["magma"],
    "pops": ["magma"],
    "flames": [
        "annot_ldblock", "formatter", "finemap",
        "ld_clump", "magma", "magmacovar", "pops"
    ],
    "heritability": ["formatter"],
    "manhattan": ["formatter"],
    "qc": [],
}

# =====================================================================
# MODULE → PARSERS
# =====================================================================

MODULE_PARSERS = {
    "annot_ldblock": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_genomeversion_parser, get_common_out_parser,
        get_annot_ldblock_parser
    ],
    "formatter": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_formatter_parser
    ],
    "finemap": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_finemap_method_parser,
        get_common_susie_arguments
    ],
    "ld_clump": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_plink_binary_parser
    ],
    "magma": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_magma_assoc_parser,
        get_magma_binary_parser
    ],
    "magmacovar": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_magma_covar_parser,
        get_magma_binary_parser
    ],
    "pops": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_pops_parser
    ],
    "flames": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_flames_common_parser,
        get_common_pops_parser, get_common_magma_assoc_parser
    ],
    "sumstat_filter": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_sumstat_filter_parser
    ],
    "imputation": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_imputation_parser
    ],
    "heritability": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_ldsc_common_parser
    ],
    "manhattan": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_assoc_plot_parser
    ],
    "qc_summary": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_qc_parser
    ],
}

# =====================================================================
# MODULE → WORKFLOW FUNCTION
# =====================================================================

WORKFLOW_FUNCS = {
    "annot_ldblock": run_annot_ldblock,
    "formatter": run_formatter_direct,
    "sumstat_filter": run_sumstat_filter_direct,
    "imputation": run_sumstat_imputation_direct,
    "finemap": run_susie_direct,
    "ld_clump": run_ld_clump_direct,
    "magma": run_magma_direct,
    "magmacovar": run_magma_covar_direct,
    "pops": run_pops_direct,
    "flames": run_flames_direct,
    "heritability": run_ldsc_direct,
    "manhattan": run_assoc_plot_direct,
    "qc_summary": run_qc_summary_direct,
}

# =====================================================================
# DAG RESOLVER
# =====================================================================

def resolve_dependencies(mods):
    resolved, visited = [], set()

    def dfs(x):
        if x in visited:
            return
        visited.add(x)
        for dep in MODULE_DEPENDENCIES.get(x, []):
            dfs(dep)
        resolved.append(x)

    for m in mods:
        dfs(m)

    return resolved


# =====================================================================
# PRINT HELP BLOCK FOR MISSING ARGUMENTS
# =====================================================================

def print_full_pipeline_help(modules, parser):
    print("\n❗ Required arguments are missing for the selected modules.")
    print("   Displaying full pipeline help:\n")

    print("📘 Detailed Pipeline Execution Plan:\n")
    for idx, m in enumerate(modules, 1):
        print(f" {idx}) {m}")
        for line in MODULE_DESCRIPTIONS.get(m, ["No description available."]):
            print(f"      • {line}")
        print("")
    print("👇 Below are the arguments relevant to this pipeline:\n")
    parser.print_help()


# =====================================================================
# MAIN
# =====================================================================

def main():

    # -------------------------------------------------------------
    # Stage 1 — Minimal parser to detect --modules and flags
    # -------------------------------------------------------------
    mini = argparse.ArgumentParser(add_help=False)
    mini.add_argument("--modules", nargs="*")
    mini.add_argument("--apply-filter", action="store_true")
    mini.add_argument("--apply-imputation", action="store_true")
    mini.add_argument("--apply-manhattan", action="store_true")
    mini.add_argument("-h", "--help", action="store_true")

    a1, _ = mini.parse_known_args()

    # -------------------------------------------------------------
    # If no modules given → show module table
    # -------------------------------------------------------------
    if not a1.modules:
        from rich.console import Console
        from rich.table import Table

        console = Console()

        console.print("""
Each module operates directly on your harmonised GWAS VCF and required
reference resources. No manual intermediate file inputs are required.
""")
        
        console.print("""
[bold cyan]Usage Example[/bold cyan]
  postgwas pipeline --modules finemap --apply-filter --apply-imputation
""")

        console.print("""
[bold cyan]PostGWAS Pipeline — Available Modules[/bold cyan] """)


        t = Table(title=" ", header_style="bold cyan", show_lines=True)
        t.add_column("Module", style="magenta", no_wrap=True)
        t.add_column("Description", style="white")

        for k in MODULE_DESCRIPTIONS:
            t.add_row(k, MODULE_DESCRIPTIONS[k][0])

        console.print(t)

        opt = Table(title="Optional Flags (Act as Modules)", show_lines=True, header_style="bold green")
        opt.add_column("Flag", style="yellow", no_wrap=True)
        opt.add_column("Description")

        opt.add_row("--apply-filter", "Add QC filtering before downstream modules.")
        opt.add_row("--apply-imputation", "Add PRED-LD imputation step.")
        opt.add_row("--apply-manhattan", "Add Manhattan/QQ plotting step.")
        opt.add_row("--heritability", "Run LDSC heritability analysis.")

        console.print(opt)
        sys.exit(0)

    # -------------------------------------------------------------
    # Expand optional flags into pipeline modules
    # -------------------------------------------------------------
    modules = list(a1.modules)

    if a1.apply_filter and "sumstat_filter" not in modules:
        modules.insert(0, "sumstat_filter")

    if a1.apply_imputation and "imputation" not in modules:
        if "sumstat_filter" in modules:
            modules.insert(modules.index("sumstat_filter") + 1, "imputation")
        else:
            modules.insert(0, "imputation")

    if a1.apply_manhattan and "manhattan" not in modules:
        modules.append("manhattan")

    # -------------------------------------------------------------
    # Resolve DAG dependencies
    # -------------------------------------------------------------
    full_modules = resolve_dependencies(modules)

    # -------------------------------------------------------------
    # Stage 2 — Build dynamic combined parser
    # -------------------------------------------------------------
    parent_parsers = []
    seen = set()

    for m in full_modules:
        for fn in MODULE_PARSERS[m]:
            if fn not in seen:
                p = fn()
                # remove built-in help
                for act in list(p._actions):
                    if "-h" in act.option_strings or "--help" in act.option_strings:
                        p._actions.remove(act)
                parent_parsers.append(p)
                seen.add(fn)

    parser = HelpOnErrorArgumentParser(
        prog="postgwas pipeline",
        description="PostGWAS pipeline",
        formatter_class=RichHelpFormatter,
        parents=parent_parsers
    )

    parser.add_argument("--modules", nargs="*")
    parser.add_argument("--apply-filter", action="store_true")
    parser.add_argument("--apply-imputation", action="store_true")
    parser.add_argument("--apply-manhattan", action="store_true")
    parser.add_argument("--heritability", action="store_true")

    # -------------------------------------------------------------
    # HELP MODE
    # -------------------------------------------------------------
    if a1.help:
        print_full_pipeline_help(full_modules, parser)
        sys.exit(0)

    # -------------------------------------------------------------
    # FULL PARSE — Catch missing required arguments
    # -------------------------------------------------------------
    try:
        args = parser.parse_args()
    except ValueError as e:
        if str(e) == "ARGPARSE_ERROR":
            print_full_pipeline_help(full_modules, parser)
            sys.exit(2)
        else:
            raise

    # -------------------------------------------------------------
    # EXECUTION
    # -------------------------------------------------------------
    print("\n🧬 Pipeline execution order:")
    for m in full_modules:
        print(f"  → {m}")

    for m in full_modules:
        print(f"\n🚀 Running module: {m}")
        WORKFLOW_FUNCS[m](args)

    print("\n🎉 Pipeline complete.\n")


if __name__ == "__main__":
    main()
