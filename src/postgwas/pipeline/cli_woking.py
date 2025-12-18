#!/usr/bin/env python3
"""
PostGWAS — Pipeline CLI (v31, Implicit Module Anchoring)

Changes in this version
-----------------------
• FIX: Solved ordering issue where implicit modules (like 'ld_clump' required by 'finemap')
       ran too early because they weren't explicitly listed in '--modules'.
• LOGIC: Added a 'tools_to_anchor' set that expands the module list to include 
         dependencies before applying the anchor logic.
"""

import argparse
import sys
import io
from contextlib import redirect_stdout, redirect_stderr
from rich_argparse import RichHelpFormatter
from rich.console import Console
from rich.table import Table
from postgwas.harmonisation.cli import run_harmonisation

# =====================================================================
# CUSTOM PARSER
# =====================================================================

class HelpOnErrorArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError("ARGPARSE_ERROR")


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
    get_common_magma_covar_parser,
    get_common_pops_parser,
    get_flames_common_parser,
    get_common_imputation_parser,
    get_ldsc_common_parser,
    get_annot_ldblock_parser,
    get_assoc_plot_parser,
    get_magma_binary_parser,
    get_bcftools_binary_parser,
    get_plink_binary_parser
)

def get_qc_parser():
    parser = argparse.ArgumentParser(add_help=False)
    g = parser.add_argument_group("QC Summary Arguments")
    g.add_argument("--qc_output_prefix")
    return parser


# =====================================================================
# MODULE DESCRIPTIONS
# =====================================================================

MODULE_DESCRIPTIONS = {
    "harmonisation": [
        "Standardises summary statistics and converts them to VCF format.",
        "(Runs after Imputation to prepare data for downstream tools)."
    ],
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
    "post_imputation_filter": [
        "SECONDARY FILTER: Runs automatically after Harmonisation.",
        "Ensures imputed variants meet QC thresholds before analysis."
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
        "Runs LD-clumping based on annot_ldblock ld regions to extract independent GWAS lead variants.",
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
    "qc_summary": [
        "Generates QC summaries."
    ],
}

# =====================================================================
# DAG DEPENDENCIES (CLEANED)
# =====================================================================

MODULE_DEPENDENCIES = {
    # ROOT NODES
    "annot_ldblock": [],
    "harmonisation": [],
    "sumstat_filter": [],
    
    # Placeholders
    "imputation": [],
    "post_imputation_filter": [],

    # INDEPENDENT TOOLS
    "formatter":     [],
    "magma":         [],
    "heritability":  [],
    "manhattan":     [],
    "qc_summary":    [],

    # INTERNAL DEPENDENCIES
    "ld_clump": ["annot_ldblock"], # Static dependency
    "finemap": ["ld_clump"],
    "magmacovar": ["magma"],
    "pops": ["magma"],
    
    # FLAMES defaults
    "flames": [
        "annot_ldblock", "finemap",
        "ld_clump", "magma", "magmacovar", "pops"
    ],
}

# =====================================================================
# MODULE → PARSERS
# =====================================================================

MODULE_PARSERS = {
    "harmonisation": [],
    
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
        get_common_susie_arguments,get_bcftools_binary_parser
    ],
    "ld_clump": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_plink_binary_parser,
        get_bcftools_binary_parser
    ],
    "magma": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_magma_assoc_parser,
        get_magma_binary_parser,get_bcftools_binary_parser
    ],
    "magmacovar": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_magma_assoc_parser,
        get_common_magma_covar_parser,
        get_magma_binary_parser,get_bcftools_binary_parser
    ],
    "pops": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_pops_parser,
        get_common_magma_assoc_parser,get_common_magma_covar_parser,
        get_bcftools_binary_parser,get_magma_binary_parser
    ],
    "flames": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, 
        get_common_magma_assoc_parser,get_common_magma_covar_parser,
        get_common_pops_parser, get_flames_common_parser,
        get_bcftools_binary_parser,get_magma_binary_parser
    ],
    "sumstat_filter": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_sumstat_filter_parser,
        get_bcftools_binary_parser
    ],
    "post_imputation_filter": [
        get_defaultresourse_parser, get_inputvcf_parser,
        get_common_out_parser, get_common_sumstat_filter_parser,
        get_bcftools_binary_parser
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

WORKFLOW_FUNCS = {
    "harmonisation": run_harmonisation,
    "annot_ldblock": run_annot_ldblock,
    "formatter": run_formatter_direct,
    "sumstat_filter": run_sumstat_filter_direct,
    "post_imputation_filter": run_sumstat_filter_direct,
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
# DAG RESOLVER & HELPERS
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

def print_full_pipeline_help(modules, parser):
    print("\n❗ Required arguments, resources, or dependencies are missing.")
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
    # Stage 1 — Minimal parser
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
        console = Console()
        console.print("\nEach module operates directly on your harmonised GWAS VCF.\n")
        console.print("[bold cyan]Usage Example[/bold cyan]\n  postgwas pipeline --modules finemap --apply-filter\n")
        
        t = Table(title="PostGWAS Pipeline — Available Modules", header_style="bold cyan", show_lines=True)
        t.add_column("Module", style="magenta", no_wrap=True)
        t.add_column("Description", style="white")

        for k in MODULE_DESCRIPTIONS:
            t.add_row(k, MODULE_DESCRIPTIONS[k][0])

        console.print(t)
        sys.exit(0)

    # -------------------------------------------------------------
    # Expand optional flags into pipeline modules
    # -------------------------------------------------------------
    modules = list(a1.modules)
    
    # 1. Handle Filter (First Pass)
    if a1.apply_filter and "sumstat_filter" not in modules:
        modules.insert(0, "sumstat_filter")
        
    # 2. Handle Imputation
    if a1.apply_imputation and "imputation" not in modules:
        if "sumstat_filter" in modules:
            modules.insert(modules.index("sumstat_filter") + 1, "imputation")
        else:
            modules.insert(0, "imputation")
            
    # 3. Handle Manhattan
    if a1.apply_manhattan and "manhattan" not in modules:
        modules.append("manhattan")

    # =========================================================
    # 🚀 CHAIN BUILDER LOGIC (Universal)
    # =========================================================
    
    # 1. Start of the Chain
    current_parent = "annot_ldblock"
    
    # 2. Filter (Pass 1)
    if "sumstat_filter" in modules:
        if current_parent not in MODULE_DEPENDENCIES["sumstat_filter"]:
            MODULE_DEPENDENCIES["sumstat_filter"].append(current_parent)
        current_parent = "sumstat_filter"
        
    # 3. Imputation Block
    if "imputation" in modules:
        if current_parent not in MODULE_DEPENDENCIES["imputation"]:
            MODULE_DEPENDENCIES["imputation"].append(current_parent)
        
        if "imputation" not in MODULE_DEPENDENCIES["harmonisation"]:
            MODULE_DEPENDENCIES["harmonisation"].append("imputation")
        if "harmonisation" not in modules:
            modules.append("harmonisation")
            
        current_parent = "harmonisation"

        if a1.apply_filter:
            if "harmonisation" not in MODULE_DEPENDENCIES["post_imputation_filter"]:
                MODULE_DEPENDENCIES["post_imputation_filter"].append("harmonisation")
            if "post_imputation_filter" not in modules:
                modules.append("post_imputation_filter")
                
            current_parent = "post_imputation_filter"

    # 4. UNIVERSAL ANCHOR (UPDATED for Implicit Modules)
    #    We identify which tools are 'active' (even if implicitly) and anchor them.
    
    # Define set of active tools to anchor
    tools_to_anchor = set()
    for m in modules:
        tools_to_anchor.add(m)
    
    # Add implicit dependencies (Known Logic)
    if "finemap" in modules:
        tools_to_anchor.add("ld_clump")
    if "magmacovar" in modules:
        tools_to_anchor.add("magma")
    if "pops" in modules:
        tools_to_anchor.add("magma")
        
    analysis_tools = [
        "annot_ldblock", 
        "heritability", 
        "manhattan", 
        "qc_summary", 
        "magma", 
        "ld_clump", 
        "finemap"
    ]
    
    for tool in analysis_tools:
        if tool in tools_to_anchor:
            # Anchor only if strictly necessary to avoid cycles
            if tool != current_parent:
                if current_parent not in MODULE_DEPENDENCIES[tool]:
                    MODULE_DEPENDENCIES[tool].append(current_parent)

    # 5. FLAMES Internal Chain
    if "flames" in modules:
        chain_pairs = [
            ("ld_clump", "finemap"),
            ("finemap", "magma"),
            ("magma", "magmacovar"),
            ("magmacovar", "pops")
        ]
        for parent, child in chain_pairs:
            if parent not in MODULE_DEPENDENCIES[child]:
                 MODULE_DEPENDENCIES[child].append(parent)

    full_modules = resolve_dependencies(modules)

    # -------------------------------------------------------------
    # Stage 2 — Build dynamic combined parser
    # -------------------------------------------------------------
    parent_parsers = []
    seen = set()

    for m in full_modules:
        for fn in MODULE_PARSERS[m]:
            if fn not in seen:
                parent_parsers.append(fn())
                seen.add(fn)

    custom_usage = "postgwas pipeline [options] --modules " + " ".join(full_modules)
    
    parser = HelpOnErrorArgumentParser(
        prog="postgwas pipeline",
        formatter_class=RichHelpFormatter,
        parents=parent_parsers,
        conflict_handler='resolve',
        usage=custom_usage
    )

    parser.add_argument("--modules", nargs="*")
    parser.add_argument("--apply-filter", action="store_true")
    parser.add_argument("--apply-imputation", action="store_true")
    parser.add_argument("--apply-manhattan", action="store_true")
    parser.add_argument("--heritability", action="store_true")

    if a1.help:
        print_full_pipeline_help(full_modules, parser)
        sys.exit(0)

    # -------------------------------------------------------------
    # FULL PARSE — Catch ARGPARSE errors
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
    # VALIDATION — Check Required Arguments
    # -------------------------------------------------------------
    REQUIRED_ARGS = {
        "annot_ldblock": ["vcf", "genome-version","sample_id","outdir","ld-region-dir"],
        "sumstat_filter": ["vcf","sample_id","outdir"],
        "post_imputation_filter": ["vcf","sample_id","outdir"], 
        "imputation": ["vcf","sample_id","outdir","ref_ld","gwas2vcf_resource"],
        "formatter": ["vcf","sample_id","outdir"],
        "finemap": ["vcf","sample_id","outdir","ld-region-dir","finemap_ld_ref"],
        "ld_clump": ["vcf", "genome-version","sample_id","outdir","ld-region-dir"],
        "magma": ["vcf","sample_id","outdir","ld_ref","gene_loc_file"],
        "magmacovar": ["vcf","sample_id","outdir","ld_ref","gene_loc_file","covariates"],
        "pops": ["vcf","sample_id","outdir","ld_ref","feature_mat_prefix","pops_gene_loc_file"],
        "flames": ["vcf","sample_id","outdir","ld_ref","covariates","feature_mat_prefix","pops_gene_loc_file","flames_annot_dir"],
        "heritability": ["vcf","sample_id","outdir","merge-alleles","ref-ld-chr","w-ld-chr"],
        "manhattan": ["vcf","sample_id","outdir"],
        "qc_summary": ["vcf","sample_id","outdir"],
    }

    missing_errors = []
    
    for m in full_modules:
        if m in REQUIRED_ARGS:
            for arg_name in REQUIRED_ARGS[m]:
                attr_name = arg_name.replace("-", "_")
                val = getattr(args, attr_name, None)
                if not val:
                    missing_errors.append(f"Module '[bold cyan]{m}[/bold cyan]' requires argument [bold red]--{arg_name}[/bold red]")

    if missing_errors:
        print("\n❌ [bold red]Missing Required Arguments:[/bold red]")
        for err in missing_errors:
            print(f"   • {err}")
        print_full_pipeline_help(full_modules, parser)
        sys.exit(2)

    # -------------------------------------------------------------
    # EXECUTION
    # -------------------------------------------------------------
    try:
        for m in full_modules:            
            f_stdout = io.StringIO()
            f_stderr = io.StringIO()
            should_print_output = True
            
            try:
                with redirect_stdout(f_stdout), redirect_stderr(f_stderr):
                    WORKFLOW_FUNCS[m](args)
            
            except (ValueError, OSError, AttributeError, TypeError, SystemExit) as inner_e:
                err_msg_inner = str(inner_e).lower()
                is_validation = False
                
                if isinstance(inner_e, SystemExit) and inner_e.code != 0:
                    is_validation = True
                elif "required" in err_msg_inner:
                    is_validation = True
                elif "executable not found" in err_msg_inner:
                    is_validation = True
                elif isinstance(inner_e, AttributeError) and "'Namespace'" in str(inner_e):
                    is_validation = True
                elif isinstance(inner_e, TypeError) and "nonetype" in err_msg_inner:
                    is_validation = True
                
                if is_validation:
                    should_print_output = False
                
                raise
            
            finally:
                if should_print_output:
                    print(f_stdout.getvalue(), end='')
                    print(f_stderr.getvalue(), end='', file=sys.stderr)

    except (ValueError, OSError, AttributeError, TypeError, SystemExit) as e:
        if isinstance(e, SystemExit):
            if e.code != 0:
                print_full_pipeline_help(full_modules, parser)
                sys.exit(e.code)
            sys.exit(0)

        err_msg = str(e).lower()
        
        if "required" in err_msg:
             print_full_pipeline_help(full_modules, parser)
             sys.exit(2)
        elif "executable not found" in err_msg or "no such file or directory" in err_msg:
             print(f"\n❌ [bold red]Error:[/bold red] {e}") 
             print_full_pipeline_help(full_modules, parser)
             sys.exit(2)
        elif isinstance(e, AttributeError) and "'Namespace' object has no attribute" in str(e):
             print_full_pipeline_help(full_modules, parser)
             sys.exit(2)
        elif isinstance(e, TypeError) and "nonetype" in err_msg:
             print_full_pipeline_help(full_modules, parser)
             sys.exit(2)

        raise

    print("\n🎉 Pipeline complete.\n")


if __name__ == "__main__":
    main()