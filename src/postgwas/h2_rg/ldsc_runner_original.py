import os
import sys
import subprocess
from pathlib import Path
import datetime

# =========================================================
# UTILS & LOGGING
# =========================================================

# ANSI Colors for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_step(msg):
    print(f"\n{Colors.HEADER}=================================================={Colors.ENDC}")
    print(f"{Colors.BOLD}STEP: {msg}{Colors.ENDC}")
    print(f"{Colors.HEADER}=================================================={Colors.ENDC}")

def print_error(step_name, specific_msg):
    print(f"\n{Colors.FAIL}❌ CRITICAL ERROR IN: {step_name}{Colors.ENDC}")
    print(f"{Colors.WARNING}Details: {specific_msg}{Colors.ENDC}")
    print(f"Check the log file for the full Docker output.")

def infer_mount_root(*paths):
    """Infer the deepest common parent directory across all input paths."""
    resolved = [str(Path(p).expanduser().resolve()) for p in paths]
    common = os.path.commonpath(resolved)
    return Path(common)

# =========================================================
# MAIN RUNNER
# =========================================================

def run_ldsc(
    sumstats_tsv: str,
    out_prefix: str,
    hm3_snplist: str,
    ldscore_dir: str,
    docker_image: str = "jibinjv/ldsc:1.0.1",
    platform: str = "linux/amd64",
    info_min: float = 0.9,
    maf_min: float = 0.01,
    samp_prev: float = 0.5,
    pop_prev: float = 0.01,
):
    # 1. Resolve Paths
    try:
        sumstats_tsv = Path(sumstats_tsv).expanduser().resolve(strict=True)
        hm3_snplist = Path(hm3_snplist).expanduser().resolve(strict=True)
        ldscore_dir = Path(ldscore_dir).expanduser().resolve(strict=True)
        out_prefix = Path(out_prefix).expanduser().resolve()
        
        # Ensure output dir exists
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
    except FileNotFoundError as e:
        print_error("PRE-FLIGHT CHECKS", f"Input file not found: {e}")
        sys.exit(1)

    # 2. Setup Logging
    log_file = out_prefix.with_suffix(".ldsc.log")
    
    def log_message(msg, to_console=False):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        formatted = f"[{timestamp}] {msg}"
        with open(log_file, "a") as f:
            f.write(formatted + "\n")
        if to_console:
            print(formatted)

    # 3. Helper to run commands with streaming output
    def run_docker_step(cmd_list, step_name):
        log_message(f"STARTING: {step_name}", to_console=False)
        log_message(f"COMMAND: {' '.join(cmd_list)}", to_console=False)
        
        print(f"{Colors.CYAN}Running {step_name}...{Colors.ENDC}")
        
        try:
            # Popen allows us to read stdout in real-time
            with subprocess.Popen(
                cmd_list, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, # Merge stderr into stdout
                text=True, 
                bufsize=1
            ) as proc:
                with open(log_file, "a") as f:
                    for line in proc.stdout:
                        print(line, end='') # Print to console
                        f.write(line)       # Write to log
                
                return_code = proc.wait()
                
                if return_code != 0:
                    raise subprocess.CalledProcessError(return_code, cmd_list)
                    
            print(f"{Colors.GREEN}✔ {step_name} Completed Successfully{Colors.ENDC}")
            
        except subprocess.CalledProcessError as e:
            log_message(f"FAILED: {step_name} (Exit Code {e.returncode})")
            print_error(step_name, f"Process exited with code {e.returncode}")
            
            # Add specific hints based on common LDSC errors
            if step_name == "MUNGE_SUMSTATS":
                print(f"\n{Colors.WARNING}HINT: Common Munge Errors:{Colors.ENDC}")
                print("1. Are your columns named correctly? (SNP, A1, A2, N, P, Z/OR/BETA)")
                print("2. Is the file comma-separated or tab-separated? (Check header)")
            
            raise e # Re-raise to stop execution

    # 4. Infer Mounts
    mount_root = infer_mount_root(sumstats_tsv, out_prefix, hm3_snplist, ldscore_dir)
    log_message(f"Mount root inferred: {mount_root}", to_console=True)

    def dpath(p: Path) -> str:
        return str(p)

    # =========================================================
    # STEP 1: MUNGE SUMSTATS
    # =========================================================
    print_step("1. Formatting Summary Statistics (Munge)")
    
    munged_prefix = out_prefix
    
    cmd_munge = [
        "docker", "run",
        f"--platform={platform}", "--rm",
        "-v", f"{mount_root}:{mount_root}",
        docker_image,
        "munge_sumstats.py",
        "--sumstats", dpath(sumstats_tsv),
        "--out", dpath(munged_prefix),
        "--merge-alleles", dpath(hm3_snplist),
        "--info-min", str(info_min),
        "--maf-min", str(maf_min)
    ]

    try:
        run_docker_step(cmd_munge, "MUNGE_SUMSTATS")
    except Exception:
        sys.exit(1) # Exit cleanly, error already printed

    # =========================================================
    # STEP 2: LDSC HERITABILITY
    # =========================================================
    print_step("2. Calculating Heritability (LDSC)")
    
    munged_file = f"{munged_prefix}.sumstats.gz"
    out_h2 = f"{munged_prefix}_Liability_scale_h2"
    
    # Verify Step 1 actually created the file
    if not Path(munged_file).exists():
        print_error("PRE-LDSC CHECK", f"Expected munged file missing: {munged_file}")
        sys.exit(1)
    
    
    # Liability-scale LDSC only if both prevalence values are provided
    if samp_prev is not None and pop_prev is not None:
        cmd_h2 = [
            "docker", "run",
            f"--platform={platform}", "--rm",
            "-v", f"{mount_root}:{mount_root}",
            docker_image,
            "ldsc.py",
            "--h2", munged_file,
            "--ref-ld-chr", dpath(ldscore_dir) + "/",
            "--w-ld-chr", dpath(ldscore_dir) + "/",
            "--out", out_h2,
            "--samp-prev", str(samp_prev),
            "--pop-prev", str(pop_prev),
        ]

        try:
            run_docker_step(cmd_h2, "LDSC_Liability_scale_HERITABILITY")
        except Exception:
            sys.exit(1)


    out_h3 = f"{munged_prefix}_h2"

    cmd_h3 = [
        "docker", "run",
        f"--platform={platform}", "--rm",
        "-v", f"{mount_root}:{mount_root}",
        docker_image,
        "ldsc.py",
        "--h2", munged_file,
        "--ref-ld-chr", dpath(ldscore_dir) + "/",
        "--w-ld-chr", dpath(ldscore_dir) + "/",
        "--out", out_h3,
    ]

    try:
        run_docker_step(cmd_h3, "LDSC_HERITABILITY")
    except Exception:
        sys.exit(1)
        
    # =========================================================
    # SUCCESS SUMMARY
    # =========================================================
    print(f"\n{Colors.GREEN}{Colors.BOLD}🎉 PIPELINE COMPLETED SUCCESSFULLY!{Colors.ENDC}")
    print(f"   • Log file:   {log_file}")
    print(f"   • Munged sum: {munged_file}")
    print(f"   • H2 Results: {out_h2}.log")