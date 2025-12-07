import subprocess
from pathlib import Path
import sys
import pandas as pd
import csv

def validate_locus_file(path):
    # --- Step 1: read first few lines to guess delimiter ---
    try:
        with open(path, "r") as f:
            sample = f.read(2048)
    except Exception as e:
        print(f"❌ ERROR: Cannot open locus file '{path}'.\n{e}")
        sys.exit(1)

    # Try Sniffer first
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", " "])
        sep = dialect.delimiter
    except Exception:
        if "\t" in sample:
            sep = "\t"
        elif "," in sample:
            sep = ","
        elif ";" in sample:
            sep = ";"
        else:
            sep = r"\s+"

    # --- Step 2: read using inferred delimiter ---
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception as e:
        print(f"❌ ERROR: Cannot parse locus file '{path}'. Using inferred separator '{sep}'.\n{e}")
        sys.exit(1)

    # --- Step 3: normalize column names ---
    df.columns = [c.strip().upper() for c in df.columns]

    # --- Step 4: validate required columns ---
    required = {"CHR", "START", "END"}
    missing = required - set(df.columns)

    if missing:
        print(
            "❌ ERROR: Locus file missing required columns.\n"
            f"   Missing: {', '.join(missing)}\n"
            f"📌 Required: CHR, START, END\n"
            f"📄 Found: {', '.join(df.columns)}\n"
            f"🔍 Inferred delimiter: '{sep}'"
        )
        sys.exit(1)

    print(f"✅ Locus file loaded successfully using delimiter '{sep}'.")
    return df



def run_susie(
    locus_file,
    sumstat_file,
    sample_id,
    ld_ref,
    plink,
    output_folder,
    lp_threshold="7.3",
    L="10",
    workers="auto",
    min_ram_per_worker_gb="4",
    timeout_ld_seconds="180",
    timeout_susie_seconds="180",
    skip_mhc=False,
    mhc_start="25000000",
    mhc_end="35000000",
    verbose=False,
):
    
    # Ensure log file parent folder exists
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    log_file = output_folder / f"{sample_id}_run_susie.log"

    script_dir = Path(__file__).resolve().parent
    rscript_file = script_dir / "run_susie_parallel_cli.R"

    if not rscript_file.exists():
        raise FileNotFoundError(f"❌ ERROR: Cannot find R script at {rscript_file}")

    cmd = [
        "Rscript", str(rscript_file),
        "--locus_file", locus_file,
        "--sumstat_file", sumstat_file,
        "--sample_id", sample_id,
        "--ld_ref", ld_ref,
        "--plink", plink,
        "--SUSIE_Analysis_folder", str(output_folder),
        "--lp_threshold", lp_threshold,
        "--L", L,
        "--workers", workers,
        "--min_ram_per_worker_gb", min_ram_per_worker_gb,
        "--timeout_ld_seconds", timeout_ld_seconds,
        "--timeout_susie_seconds", timeout_susie_seconds,
        "--mhc_start", mhc_start,
        "--mhc_end", mhc_end,
    ]

    if skip_mhc:
        cmd.append("--skip_mhc")

    if verbose:
        cmd.append("--verbose")

    # Convert all items to string
    cmd = list(map(str, cmd))

    # Redirect stdout + stderr to log file
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)
