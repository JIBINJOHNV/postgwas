"""
PostGWAS Pipeline Runners (v56, Dynamic Directories & Logic Fixes)

Changes:
• DYNAMIC DIRS: Directories now use step numbers from execution.py (e.g., 01_formatter, 06_formatter).
• MAGMA: Auto-calculates batches and runs annotation first.
• CONTEXT: Fixes dictionary key errors in Magma/Covar/Finemap/PoPS.
"""

import os
import shutil
import sys
import traceback
import multiprocessing
import pandas as pd  # Required for MAGMA batch calculation

# Workflow Imports
from postgwas.annot_ldblock.workflows import run_annot_ldblock
from postgwas.formatter.workflows import run_formatter_direct
from postgwas.sumstat_filter.workflows import run_sumstat_filter_direct
from postgwas.imputation.workflows import run_sumstat_imputation_direct
from postgwas.finemap.workflows import run_susie_direct
from postgwas.ld_clump.workflows import run_ld_clump_direct
# Note: Added run_magma_annot_direct
from postgwas.gene_assoc.workflows import run_magma_direct
from postgwas.magmacovar.workflows import run_magma_covar_direct
from postgwas.pops.workflows import run_pops_direct
from postgwas.flames.workflows import run_flames_direct
from postgwas.h2_rg.workflows import run_ldsc_direct
from postgwas.manhattan.workflows import run_assoc_plot_direct
from postgwas.qc_summary.workflows import run_qc_summary_direct

# ============================================================
# HELPER: Directory Manager (DYNAMIC)
# ============================================================

def setup_subdir(args, base_name):
    """
    Dynamically creates subdirectories based on execution order.
    Example: If base_name is 'formatter' and step is 2, creates '02_formatter'.
    """
    # 1. Get step number from args (injected by execution.py)
    # Default to '00' if running standalone or if _step_num is missing
    prefix = getattr(args, "_step_num", "00")
    
    # 2. Construct dynamic folder name
    folder_name = f"{prefix}_{base_name}"
    
    original_outdir = args.outdir
    new_outdir = os.path.join(original_outdir, folder_name)
    os.makedirs(new_outdir, exist_ok=True)
    
    # Update args.outdir temporarily for the runner
    args.outdir = new_outdir
    return original_outdir

# ============================================================
# HELPER: Batch Calculator (MAGMA)
# ============================================================

def get_optimal_magma_batches(annot_file, user_requested_batches=None):
    """
    Calculates the optimal number of MAGMA batches based on Data and CPU limits.
    Prevents 'batch X is empty' errors on small datasets.
    """
    sys_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
    
    try:
        # Read only the chromosome column (index 1) to verify unique count
        df = pd.read_csv(annot_file, delim_whitespace=True, header=None, usecols=[1])
        data_chroms = df[1].nunique()
    except Exception:
        # Fallback if file read fails
        data_chroms = 22 
        
    user_limit = int(user_requested_batches) if user_requested_batches else 999
    
    # Logic: Batches cannot exceed chromosome count
    optimal = min(data_chroms, sys_cpus, user_limit)
    optimal = max(1, optimal) # Ensure at least 1

    print(f"   ℹ️  [MAGMA] Auto-adjusted batches to: {optimal} (Chromosomes: {data_chroms})")
    return optimal

# ============================================================
# HELPER: Dependency Debugger
# ============================================================

def check_and_resolve_binaries(args, required_tools):    
    def inject_to_path(binary_path):
        directory = os.path.dirname(binary_path)
        if directory not in os.environ["PATH"]:
            os.environ["PATH"] = directory + os.pathsep + os.environ["PATH"]

    if "bcftools" in required_tools or hasattr(args, "bcftools"):
        cmd = getattr(args, 'bcftools', 'bcftools')
        resolved = shutil.which(cmd)
        if resolved:
            args.bcftools = resolved
            inject_to_path(resolved)
        else:
            raise FileNotFoundError(f"Missing executable: {cmd}")

    if "plink" in required_tools:
        cmd = getattr(args, 'plink', 'plink')
        resolved = shutil.which(cmd)
        if not resolved and cmd == "plink":
             resolved = shutil.which("plink2")
        
        if resolved:
            args.plink = resolved
            inject_to_path(resolved)
        else:
            raise FileNotFoundError(f"Missing executable: {cmd}")

    for tool in ["tabix", "bgzip"]:
        if tool in required_tools:
            path = shutil.which(tool)
            if path:
                inject_to_path(path)
            else:
                raise FileNotFoundError(f"Missing tool: {tool}")

# ============================================================
# HELPER: VCF Extractor
# ============================================================

def get_vcf_from_context(ctx):
    if "post_imputation_filter" in ctx:
        return ctx["post_imputation_filter"]["filtered_vcf"]
    
    if "imputation" in ctx:
        data = ctx["imputation"]
        if isinstance(data, dict):
            return data.get("GRCh38", data.get("GRCh37"))
        return data 

    if "sumstat_filter" in ctx:
        return ctx["sumstat_filter"]["filtered_vcf"]

    if "annot_ldblock" in ctx:
        return ctx["annot_ldblock"]["annotated_vcf"]
        
    return None

# ============================================================
# RUNNERS
# ============================================================

def run_annot_ldblock_runner(args, ctx):
    root = setup_subdir(args, "annot_ldblock")
    # vcf = get_vcf_from_context(ctx)
    # if vcf: args.vcf = vcf
    try:
        outputs = run_annot_ldblock(args)
        ctx["annot_ldblock"] = outputs
    finally:
        args.outdir = root
        args.vcf = outputs["annotated_vcf"]
    return outputs

def run_sumstat_filter_runner(args, ctx):
    # Determine base name based on context
    if "imputation" in ctx:
        base_name = "filter_post_imp"
    else:
        base_name = "filter_pre_imp"

    root = setup_subdir(args, base_name)
    # vcf = get_vcf_from_context(ctx)
    # if vcf: args.vcf = vcf
    try:
        outputs = run_sumstat_filter_direct(args)
        if "imputation" in ctx:
            ctx["post_imputation_filter"] = outputs
        else:
            ctx["sumstat_filter"] = outputs
    finally:
        args.outdir = root
        args.vcf = outputs["filtered_vcf"]
    return outputs

def run_formatter_runner(args, ctx):
    root = setup_subdir(args, "formatter")

    # 1. Determine Formats
    if not hasattr(args, "format") or args.format is None:
        active_formats = set()
    else:
        active_formats = set(args.format)

    active_modules = set(args.modules) if args.modules else set()

    if args.apply_imputation or "imputation" in active_modules:
        active_formats.add("ldpred")
    if "heritability" in active_modules:
        active_formats.add("ldsc")
    if any(m in active_modules for m in ["magma", "magmacovar", "pops", "flames"]):
        active_formats.add("magma")
    if "finemap" in active_modules or "flames" in active_modules:
        active_formats.add("finemap")

    args.format = list(active_formats)

    # 2. Resolve Binaries
    required = ["bcftools", "tabix", "bgzip"]
    if any(fmt in args.format for fmt in ["ldpred", "magma", "finemap"]):
        required.append("plink")
        
    check_and_resolve_binaries(args, required)

    # 3. Determine Input VCF
    # vcf = get_vcf_from_context(ctx)
    # if vcf: args.vcf = vcf

    if args.vcf and not os.path.exists(args.vcf):
        print(f"\n❌ [bold red]File Not Found:[/bold red] The input VCF does not exist:\n   {args.vcf}")
        raise FileNotFoundError(f"Input VCF missing: {args.vcf}")

    try:
        outputs = run_formatter_direct(args)
        ctx["formatter"] = outputs
    except Exception:
        print("\n\n🔥 [CRITICAL FORMATTER FAILURE] 🔥")
        traceback.print_exc()
        print("---------------------------------------------------")
        raise
    finally:
        args.outdir = root
    return outputs

def run_imputation_runner(args, ctx):
    root = setup_subdir(args, "imputation")
    
    if "ldpred" in ctx["formatter"] and ctx["formatter"]["ldpred"]:
        ldpred_folder = ctx["formatter"]["ldpred"].get("ldpred_folder")
        if not ldpred_folder:
             raise ValueError("Formatter did not return a valid 'ldpred_folder' path.")
    else:
        raise ValueError("Imputation requires 'ldpred' format from formatter.")

    args.predld_input_dir = ldpred_folder
    
    try:
        outputs = run_sumstat_imputation_direct(args)
        args.vcf = outputs["GRCh37"]
    finally:
        args.outdir = root
    return outputs

def run_finemap_runner(args, ctx):
    root = setup_subdir(args, "finemap")
    args.finemap_input_file = ctx["formatter"]["finemap"]["susie_input"]
    args.locus_file = ctx["ld_clump"]["ldpruned_sig_file"]
    try:
        outputs = run_susie_direct(args)
        print("finemap returned:")
        print(outputs)
        # FIX: Save entire output dict
        ctx["finemap"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_ld_clump_runner(args, ctx):
    root = setup_subdir(args, "ld_clump")
    #vcf = get_vcf_from_context(ctx)
    args.ld_mode = "by_regions"
    # if vcf: args.vcf = vcf
    try:
        outputs = run_ld_clump_direct(args)
        ctx["ld_clump"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_magma_runner(args, ctx):
    root = setup_subdir(args, "magma")
    
    magma_inputs = ctx["formatter"]["magma"]
    args.snp_loc_file = magma_inputs["snp_loc_file"]
    args.pval_file = magma_inputs["pval_file"]

    try:
        # 1. RUN ANALYSIS
        outputs = run_magma_direct(args, ctx)
        ctx["magma"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_magmacovar_runner(args, ctx):
    root = setup_subdir(args, "magma_covar")
    # FIX: Correct key is "magma", not "magma_gene"
    args.magama_gene_assoc_raw = ctx["magma"]["magma_genes_raw"]
    try:
        outputs = run_magma_covar_direct(args)
        ctx["magma_covar"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_pops_runner(args, ctx):
    root = setup_subdir(args, "pops")
    
    if "magma" in ctx and "magma_genes_prefix" in ctx["magma"]:
        args.magma_assoc_prefix = ctx["magma"]["magma_genes_prefix"]
    else:
        print("⚠️ Warning: MAGMA prefix not found in context. PoPS might fail.")
        
    args.pops_verbose = True
    
    try:
        # FIX: Pass 'ctx' so workflows can write to it
        outputs = run_pops_direct(args, ctx)
        
        if "pops_output" not in ctx and isinstance(outputs, dict):
             ctx["pops_output"] = outputs.get('pops_file')
    finally:
        args.outdir = root
    return outputs

def run_flames_runner(args, ctx):
    root = setup_subdir(args, "flames")
    # FIX: Ensure keys exist in context
    args.finemap_cred_dir = ctx["finemap"]["flames_input"]
    args.magma_genes_out = ctx["magma"]["magma_genes_out"]
    args.magma_tissue_covar_results = ctx["magma_covar"]
    args.pops_score_file = ctx["pops_output"]
    try:
        outputs = run_flames_direct(args)
        ctx["flames"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_heritability_runner(args, ctx):
    root = setup_subdir(args, "heritability")
    args.ldsc_inut = ctx["formatter"]["ldsc"]["ldsc_file"]
    try:
        # FIX: Pass file explicitly to updated run_ldsc_direct
        outputs = run_ldsc_direct(args)
        ctx["heritability"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_manhattan_runner(args, ctx):
    root = setup_subdir(args, "manhattan")
    # vcf = get_vcf_from_context(ctx)
    # if vcf: args.vcf = vcf
    try:
        outputs = run_assoc_plot_direct(args)
        ctx["manhattan"] = outputs
    finally:
        args.outdir = root
    return outputs

def run_qc_summary_runner(args, ctx):
    root = setup_subdir(args, "qc_summary")
    # vcf = get_vcf_from_context(ctx)
    # if vcf: args.vcf = vcf
    try:
        outputs = run_qc_summary_direct(args)
        ctx["qc_summary"] = outputs
    finally:
        args.outdir = root
    return outputs

# ============================================================
# PIPELINE REGISTRY
# ============================================================

RUNNERS = {
    "annot_ldblock": run_annot_ldblock_runner,
    "sumstat_filter": run_sumstat_filter_runner,
    "post_imputation_filter": run_sumstat_filter_runner, 
    "formatter": run_formatter_runner,
    "imputation": run_imputation_runner,
    "finemap": run_finemap_runner,
    "ld_clump": run_ld_clump_runner,
    "magma": run_magma_runner,
    "magmacovar": run_magmacovar_runner,
    "pops": run_pops_runner,
    "flames": run_flames_runner,
    "heritability": run_heritability_runner,
    "manhattan": run_manhattan_runner,
    "qc_summary": run_qc_summary_runner,
}