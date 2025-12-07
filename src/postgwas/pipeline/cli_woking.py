#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v9)

New Feature
-----------
When user runs:
    postgwas pipeline --modules X [...] --help

A detailed pipeline execution plan is shown BEFORE the argument help:

    1) annot_ldblock
          • Annotates variants with LD blocks...
          • Requires: ...

    2) formatter
          • Converts harmonised VCF...

    ...

This is fully dynamic and automatically includes:
    --apply-filter
    --apply-imputation
    --apply-manhattan
"""

import argparse
import sys
from rich_argparse import RichHelpFormatter

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
    group = parser.add_argument_group("QC Summary Arguments")
    group.add_argument("--qc_output_prefix")
    return parser

# =====================================================================
# MODULE DESCRIPTIONS (Dynamic Help Preview)
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
        "Requires: PLINK binary files and formatted association statistics."
    ],
    "magma": [
        "Runs MAGMA gene-level association analysis.",
        "Requires: MAGMA binary, gene annotation, formatted inputs."
    ],
    "magmacovar": [
        "Runs MAGMA gene-property (covariate) analysis.",
        "Requires: .genes.raw file and covariate matrix file."
    ],
    "pops": [
        "Computes PoPS gene prioritisation scores.",
        "Requires: MAGMA gene-level outputs + PoPS resource files."
    ],
    "flames": [
        "Runs FLAMES integrative effector-gene scoring.",
        "Combines fine-mapping, MAGMA, PoPS, annotations for final scoring."
    ],
    "heritability": [
        "Runs LDSC heritability and genetic correlation estimation.",
        "Requires: LDSC reference LD scores & weight files."
    ],
    "manhattan": [
        "Generates Manhattan & QQ association plots.",
        "Requires: formatted summary statistics file."
    ],
    "qc": [
        "Generates QC summary tables for harmonisation/filtering steps.",
        "Includes variant counts, ts/tv, SNP/INDEL ratios, allele mismatches."
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
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_genomeversion_parser,
        get_common_out_parser,
        get_annot_ldblock_parser,
    ],
    "formatter": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_formatter_parser,
    ],
    "finemap": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_finemap_method_parser,
        get_common_susie_arguments,
    ],
    "ld_clump": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_plink_binary_parser,
    ],
    "magma": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_common_magma_assoc_parser,
        get_magma_binary_parser,
    ],
    "magmacovar": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_common_magma_covar_parser,
        get_magma_binary_parser,
    ],
    "pops": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_common_pops_parser,
    ],
    "flames": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_flames_common_parser,
        get_common_pops_parser,
        get_common_magma_assoc_parser,
    ],
    "sumstat_filter": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_common_sumstat_filter_parser,
    ],
    "imputation": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_common_imputation_parser,
    ],
    "heritability": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_ldsc_common_parser,
    ],
    "manhattan": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_assoc_plot_parser,
    ],
    "qc_summary": [
        get_defaultresourse_parser,
        get_inputvcf_parser,
        get_common_out_parser,
        get_qc_parser,
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
    resolved = []
    visited = set()

    def dfs(m):
        if m in visited:
            return
        visited.add(m)
        for dep in MODULE_DEPENDENCIES.get(m, []):
            dfs(dep)
        resolved.append(m)

    for m in mods:
        dfs(m)
    return resolved

# =====================================================================
# MAIN
# =====================================================================
def main():

    # ---------------- Stage 1: Pre-parse ----------------
    s1 = argparse.ArgumentParser(add_help=False)
    s1.add_argument("--modules", nargs="*")
    s1.add_argument("--apply-filter", action="store_true")
    s1.add_argument("--apply-imputation", action="store_true")
    s1.add_argument("--apply-manhattan", action="store_true")
    s1.add_argument("-h", "--help", action="store_true")

    a1, unknown = s1.parse_known_args()

    # =================================================================
    # NEW MODULE LISTING TABLE (OPTION A) — ONLY BLOCK THAT CHANGED
    # =================================================================
    if not a1.modules:
        from rich.console import Console
        from rich.table import Table

        console = Console()
        table = Table(title="PostGWAS Pipeline Modules", header_style="bold cyan", show_lines=True)
        console.print("""
[bold cyan]PostGWAS Pipeline — Available Modules[/bold cyan]

Each module performs a specific post-GWAS analysis step using your
harmonised summary statistics VCF (plus module-specific reference files).  
All intermediate inputs (MAGMA tables, FINEMAP inputs, LDSC files, clumped sets)
are generated internally — you never provide them manually.

Below is a list of all modules and what each one does:
        """)
        table.add_column("Module", style="magenta", no_wrap=True)
        table.add_column("Description", style="white")

        table.add_row("annot_ldblock", "Annotate variants with LD-block identifiers.")
        table.add_row("formatter", "Convert harmonised VCF to tool-specific formats.")
        table.add_row("sumstat_filter", "Summary statistics QC.")
        table.add_row("imputation", "PredLD summary-statistics imputation.")
        table.add_row("finemap", "SuSiE fine-mapping.")
        table.add_row("ld_clump", "PLINK LD-clumping.")
        table.add_row("magma", "MAGMA gene-level association.")
        table.add_row("magmacovar", "MAGMA gene-property (covariate) analysis.")
        table.add_row("pops", "PoPS gene prioritisation.")
        table.add_row("flames", "Integrative effector-gene scoring.")
        table.add_row("heritability", "LDSC heritability / genetic correlation.")
        table.add_row("manhattan", "Manhattan & QQ plots.")
        table.add_row("qc_summary", "Harmonisation/filtering QC summaries.")

        console.print(table)

        opt = Table(title="Optional Flags (behave as modules)", header_style="bold green", show_lines=True)
        opt.add_column("Flag", style="yellow", no_wrap=True)
        opt.add_column("Action")

        opt.add_row("--apply-filter", "Add summary-statistics QC step.")
        opt.add_row("--apply-imputation", "Add summary-statistics imputation.")
        opt.add_row("--apply-manhattan", "Add Manhattan/QQ plotting step.")
        opt.add_row("--heritability", "Run LDSC heritability.")

        console.print(opt)

    # ---------------------- USAGE EXAMPLES ------------------------
        console.print("""
        [bold cyan]Usage Examples[/bold cyan]

        • Run LD-block annotation + formatter:
            [yellow]postgwas pipeline --modules annot_ldblock formatter[/yellow]

        • Run full gene-based pipeline (MAGMA):
            [yellow]postgwas pipeline --modules annot_ldblock formatter magma[/yellow]

        • Run fine-mapping with filter + imputation:
            [yellow]postgwas pipeline --modules finemap --apply-filter --apply-imputation[/yellow]

        • Run PoPS + MAGMA + Manhattan plot:
            [yellow]postgwas pipeline --modules magma pops --apply-manhattan[/yellow]

        • Run LDSC heritability alone:
            [yellow]postgwas pipeline --modules heritability[/yellow]

        • Run a full integrative workflow:
            [yellow]postgwas pipeline --modules annot_ldblock formatter magma finemap pops flames \\
                --apply-filter --apply-imputation --apply-manhattan --heritability[/yellow]

        """)
        sys.exit(0)

    # =================================================================
    # OPTIONAL ADD-ON MODULES
    # =================================================================
    base_modules = list(a1.modules)

    if a1.apply_filter and "sumstat_filter" not in base_modules:
        base_modules.insert(0, "sumstat_filter")

    if a1.apply_imputation and "imputation" not in base_modules:
        if "sumstat_filter" in base_modules:
            idx = base_modules.index("sumstat_filter") + 1
            base_modules.insert(idx, "imputation")
        else:
            base_modules.insert(0, "imputation")

    if a1.apply_manhattan and "manhattan" not in base_modules:
        base_modules.append("manhattan")

    # DAG expand
    modules = resolve_dependencies(base_modules)

    # ---------------- Stage 2: Build dynamic parser ----------------
    parent_parsers = []
    added = set()

    for mod in modules:
        for fn in MODULE_PARSERS[mod]:
            if fn not in added:
                p = fn()
                # remove help to avoid conflicts
                for act in list(p._actions):
                    if "-h" in act.option_strings or "--help" in act.option_strings:
                        p._actions.remove(act)
                parent_parsers.append(p)
                added.add(fn)

    final = argparse.ArgumentParser(
        prog="postgwas pipeline",
        description="PostGWAS unified pipeline runner",
        formatter_class=RichHelpFormatter,
        parents=parent_parsers,
    )

    final.add_argument("--modules", nargs="*")
    final.add_argument("--apply-filter", action="store_true")
    final.add_argument("--apply-imputation", action="store_true")
    final.add_argument("--apply-manhattan", action="store_true")
    final.add_argument("--heritability", action="store_true")

    # ---------------- HELP Mode: Print Detailed Pipeline Plan ----------------
    if a1.help:
        print("\n📘 Detailed Pipeline Execution Plan (after dependency expansion):\n")

        for idx, m in enumerate(modules, start=1):
            print(f" {idx}) {m}")

            desc = MODULE_DESCRIPTIONS.get(m, ["No description available."])
            for line in desc:
                print(f"      • {line}")
            print("")

        print("👇 Below are the arguments relevant to the modules above:\n")
        final.print_help()
        sys.exit(0)

    # ---------------- Full parse ----------------
    args = final.parse_args()

    # ---------------- EXECUTION ----------------
    print("\n🧬 Pipeline execution order:")
    for m in modules:
        print(f"  → {m}")

    for m in modules:
        print(f"\n🚀 Running module: {m}")
        WORKFLOW_FUNCS[m](args)

    print("\n🎉 Pipeline complete.\n")


if __name__ == "__main__":
    main()
