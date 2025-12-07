#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from datetime import datetime
import subprocess
from rich_argparse import RichHelpFormatter





# =====================================================================
#   RUN MODES
# =====================================================================
def run_assoc_plot_direct(args):
    """Executes assoc_plot.R in direct mode."""
    script_dir = Path(__file__).resolve().parent
    R_SCRIPT = script_dir / "assoc_plot.R"

    if not R_SCRIPT.exists():
        sys.exit(f"❌ assoc_plot.R not found in {script_dir}")

    # Build command
    cmd = ["Rscript", str(R_SCRIPT)]

    def add(flag, value):
        if value is not None:
            cmd.append(f"{flag}={value}")

    add("--genome", args.genome)
    add("--cytoband", args.cytoband)
    add("--vcf", args.vcf)
    add("--pheno", args.pheno)

    add("--nauto", args.nauto)
    add("--min-af", args.min_af)
    add("--min-lp", args.min_lp)
    add("--loglog-pval", args.loglog_pval)
    add("--cyto-ratio", args.cyto_ratio)
    add("--max-height", args.max_height)
    add("--spacing", args.spacing)

    add("--pdf", args.pdf)
    add("--png", args.png)
    add("--width", args.width)
    add("--height", args.height)
    add("--fontsize", args.fontsize)

    if args.allelic_shift:
        cmd.append("--as")
    if args.csq:
        cmd.append("--csq")

    # Determine log file
    if args.pdf:
        base = Path(args.pdf)
    elif args.png:
        base = Path(args.png)
    else:
        sys.exit("❌ Either --pdf or --png is required.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pheno_tag = args.pheno or "NO_PHENO"

    log_file = base.parent / f"{base.stem}.{pheno_tag}.assocplot_{timestamp}.log"

    # Run
    with open(log_file, "w") as log:
        log.write("COMMAND:\n" + " ".join(cmd) + "\n\n")
        log.write("------ Rscript output ------\n\n")

        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=log,
            text=True
        )

    if process.returncode != 0:
        print(f"❌ assoc_plot failed. See: {log_file}")
        sys.exit(process.returncode)

