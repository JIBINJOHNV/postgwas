import sys
from postgwas.pipeline.runners import RUNNERS

import sys
from postgwas.pipeline.runners import RUNNERS

def execute_pipeline(args, modules):
    """
    Orchestrates the execution order and manages the data context (ctx).
    Calls Runner functions from runners.py.
    """
    
    # 1. Initialize Context (Shared Memory)
    ctx = {} 
    
    execution_chain = []
    
    # Helper to add modules only if they haven't been added yet
    def add_unique(module_name):
        if module_name in RUNNERS and module_name not in execution_chain:
            execution_chain.append(module_name)

    # Helper to check if imputation is active in this run
    is_imputing = args.apply_imputation or "imputation" in modules

    # =========================================================
    # 1. PRE-PROCESSING & IMPUTATION CHAIN
    # =========================================================
    
    # A. Filter Pass 1 (Pre-Imputation)
    if args.apply_filter or "sumstat_filter" in modules:
        add_unique("sumstat_filter")

    # B. Imputation Block
    if is_imputing:
        add_unique("formatter")  # Run 1: Prepares inputs for imputation
        add_unique("imputation") # Performs Imputation
        
        # C. Filter Pass 2 (Post-Imputation)
        if args.apply_filter or "post_imputation_filter" in modules:
            add_unique("post_imputation_filter")

    # =========================================================
    # 2. ANNOTATION & FORMATTING
    # =========================================================
    
    # D. Annotation (Run AFTER Imputation so new variants get annotated)
    
    if "annot_ldblock" in modules:
        add_unique("annot_ldblock")

    # E. Formatter (Run 2: Prepares inputs for Analysis tools)
    # Logic: If we ran imputation, we usually need to re-run formatter for downstream tools.
    # We use .append() explicitly here to allow a duplicate 'formatter' step if needed.
    downstream_tools = ["finemap", "magma", "magmacovar", "pops", "flames", "ld_clump", "heritability"]
    
    if is_imputing and any(m in modules for m in downstream_tools):
        if "formatter" in execution_chain: 
            execution_chain.append("formatter") 
        else:
            add_unique("formatter")
            
    # If NOT imputing, we still need to run formatter once if analysis tools are present
    elif not is_imputing:
        if "formatter" in modules:
            add_unique("formatter")
        elif any(m in modules for m in downstream_tools):
            add_unique("formatter")

    # =========================================================
    # 3. ANALYSIS CHAIN
    # =========================================================

    if "flames" in modules:
        # FLAMES Suite dependencies
        add_unique("ld_clump")
        add_unique("finemap")
        add_unique("magma")
        add_unique("magmacovar")
        add_unique("pops")
        add_unique("flames")
    
    else:
        # Independent Tools
        if "ld_clump" in modules: 
            add_unique("ld_clump")

        if "finemap" in modules:
            add_unique("ld_clump") # Finemap needs clumping
            add_unique("finemap")

        if "magma" in modules:
            add_unique("magma")
            
        if "magmacovar" in modules:
            add_unique("magma")
            add_unique("magmacovar")

        if "pops" in modules:
            add_unique("magma")
            add_unique("magmacovar")
            add_unique("pops")

    # =========================================================
    # 4. REPORTING CHAIN
    # =========================================================
    if "heritability" in modules:
        add_unique("heritability")

    if "manhattan" in modules or args.apply_manhattan:
        add_unique("manhattan")

    if "qc_summary" in modules:
        add_unique("qc_summary")

    # =========================================================
    # 5. EXECUTION LOOP
    # =========================================================
    
    print("\n🚀 [bold green]Starting Execution Chain[/bold green]")
    
    if not execution_chain:
        print("⚠️  [bold yellow]Warning:[/bold yellow] No modules were selected for execution.")
        print("   Check that your --modules list or flags match the available steps.")

    for i, module_name in enumerate(execution_chain, 1):
        if module_name not in RUNNERS:
            continue
            
        runner_func = RUNNERS[module_name]
        print(f"   {i}. Running: [cyan]{module_name}[/cyan]")
        
        # --- NEW: Inject Step Number for Dynamic Folders ---
        # This allows runners.py to create folders like "01_sumstat_filter", "02_imputation"
        args._step_num = f"{i:02d}"
        # ---------------------------------------------------
        
        try:
            runner_func(args, ctx)
            
        except Exception as e:
            print(f"\n❌ [bold red]Pipeline Failed at step: {module_name}[/bold red]")
            raise e 
            
    print("\n✅ [bold green]All tasks completed successfully.[/bold green]\n")