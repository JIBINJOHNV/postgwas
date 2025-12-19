import argparse
import sys
from pathlib import Path
from rich_argparse import RichHelpFormatter
from postgwas.utils.main import validate_path
import shutil

# SuSiE backend (Python) ---------------------
from postgwas.finemap.main import (
    validate_locus_file,
    run_susie,
)






def run_susie_direct(args):
    print("\n       === Running Finemapping (SuSiE)  ===")

    final_plink_path = None

    # 1. Check if user provided a specific path
    if args.plink:
        if Path(args.plink).exists():
            final_plink_path = args.plink
        else:
            print(f"   ⚠️  [WARNING] The provided path does not exist: {args.plink}")
            print("       Attempting to auto-detect 'plink' in system PATH instead...")

    # 2. If we don't have a valid path yet, search PATH
    if not final_plink_path:
        # Try 'plink' first, then 'plink2'
        detected_plink = shutil.which("plink") or shutil.which("plink2")
        if detected_plink:
            print(f"   ✅  Found 'plink' in system PATH: {detected_plink}")
            final_plink_path = detected_plink

    # 3. Final Validation: If still missing, CRASH.
    if not final_plink_path:
        print("\n❌ [ERROR] PLINK executable not found!")
        print("   The pipeline failed to find PLINK in the provided argument or the system $PATH.")
        print("   Please either:")
        print("   1. Add 'plink' to your system environment, OR")
        print("   2. Provide a valid full path using --plink /path/to/plink")
        sys.exit(1)

    # FIX: Update args.plink (not args.magma)
    args.plink = final_plink_path
    
    # Validate locus file
    validate_locus_file(args.locus_file)
    
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    output_folder = Path(args.outdir)

    susie_results = run_susie(
        locus_file=args.locus_file,
        sumstat_file=args.finemap_input_file,
        sample_id=args.sample_id,
        ld_ref=args.finemap_ld_ref,
        plink=args.plink, # Now guaranteed to be valid
        output_folder=str(output_folder),
        lp_threshold=args.lp_threshold,
        L=args.L,
        workers=args.nthreads,
        min_ram_per_worker_gb=args.min_ram_per_worker_gb,
        timeout_ld_seconds=args.timeout_ld_seconds,
        timeout_susie_seconds=args.timeout_susie_seconds,
        skip_mhc=args.finemap_skip_mhc,
        mhc_start=args.finemap_mhc_start,
        mhc_end=args.finemap_mhc_end,
    )
    print(f"        👉 Results are saved in: {output_folder.resolve()}\n")
    print(f"      🎉 Finemapping using SuSiE is Completed")
    print()
    return susie_results