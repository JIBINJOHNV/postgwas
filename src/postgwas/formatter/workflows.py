#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter

# ---- Formatter engines ----
from postgwas.formatter.main import (
    create_magma_inputs,
    create_finemap_inputs,
    create_ldpred_inputs,
    create_ldsc_inputs
)



# ============================================================
#  DIRECT MODE ENGINE
# ============================================================
def run_formatter_direct(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # args.format is a LIST when nargs="+"
    selected_formats = args.format  
    # A map of format → function
    format_map = {
        "magma": create_magma_inputs,
        "finemap": create_finemap_inputs,
        "ldpred": create_ldpred_inputs,
        "ldsc": create_ldsc_inputs,
    }

    try:
        for fmt in selected_formats:
            if fmt not in format_map:
                raise ValueError(f"Unknown format: {fmt}")

            #print(f"\n➡️  Converting to {fmt.upper()} format...")
            format_map[fmt](args)
            #print(f"✅ Finished {fmt.upper()} conversion.")

    except Exception as e:
        print(f"\n❌ Formatter Failed: {e}", file=sys.stderr)
        sys.exit(1)


# ============================================================
#  PIPELINE MODE ENGINE
# ============================================================
def run_formatter_pipeline(args):
    print("\n🚀 Running Full Pipeline: Harmonisation → LD Blocks → QC Filter → Formatter")

    # ---------------------------------------------------------
    # Step 1 + Step 2 (Harmonisation + LD block annotation)
    # ---------------------------------------------------------
    #run_annot_ldblock_pipeline(args)

    # ---------------------------------------------------------
    # Step 3 (Summary Statistics QC)
    # ---------------------------------------------------------
    #run_sumstat_qc_pipeline(args)

    # ---------------------------------------------------------
    # Step 4 (Formatter)
    # ---------------------------------------------------------
    print("\n=== Running Formatter (Pipeline Mode) ===")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.format == "magma":
        create_magma_inputs(args)

    elif args.format == "finemap":
        create_finemap_inputs(args)

