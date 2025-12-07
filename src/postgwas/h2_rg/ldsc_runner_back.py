import os
import subprocess
from pathlib import Path

def infer_mount_root(*paths):
    """
    Infer the deepest common parent directory across all input paths.
    """
    resolved = [str(Path(p).expanduser().resolve()) for p in paths]
    common = os.path.commonpath(resolved)
    return Path(common)

import subprocess
from pathlib import Path
import datetime


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

    # ---------------------------------------------------------
    # Resolve absolute paths
    # ---------------------------------------------------------
    sumstats_tsv = Path(sumstats_tsv).expanduser().resolve()
    out_prefix = Path(out_prefix).expanduser().resolve()
    hm3_snplist = Path(hm3_snplist).expanduser().resolve()
    ldscore_dir = Path(ldscore_dir).expanduser().resolve()

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # Simple log file name
    # ---------------------------------------------------------
    log_file = out_prefix.with_suffix(".ldsc.log")

    def log(msg):
        with open(log_file, "a") as f:
            f.write(msg + "\n")

    # ---------------------------------------------------------
    # Infer mount root
    # ---------------------------------------------------------
    mount_root = infer_mount_root(sumstats_tsv, out_prefix, hm3_snplist, ldscore_dir)
    log(f"Mount root: {mount_root}")

    def dpath(p: Path) -> str:
        return str(p)

    def run_cmd(cmd, label):
        """Run a command and write stdout/stderr to log."""
        log(f"\n===== {label} =====")
        log("Command: " + " ".join(cmd))

        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            log("--- STDOUT ---")
            log(result.stdout)
            log("--- STDERR ---")
            log(result.stderr)
        except subprocess.CalledProcessError as e:
            log("❌ ERROR:")
            log(e.stdout or "")
            log(e.stderr or "")
            raise

    try:
        # ---------------------------------------------------------
        # Step 1 — munge_sumstats.py
        # ---------------------------------------------------------
        munged_prefix = out_prefix

        cmd_munge = [
            "docker", "run",
            f"--platform={platform}", "--rm",
            "-v", f"{mount_root}:{mount_root}",
            docker_image,
            "munge_sumstats.py",
            "--info-min", str(info_min),
            "--maf-min", str(maf_min),
            "--sumstats", dpath(sumstats_tsv),
            "--out", dpath(munged_prefix),
            "--merge-alleles", dpath(hm3_snplist),
        ]

        run_cmd(cmd_munge, "MUNGE_SUMSTATS")

        # ---------------------------------------------------------
        # Step 2 — ldsc.py --h2
        # ---------------------------------------------------------
        munged_file = f"{munged_prefix}.sumstats.gz"
        out_h2 = f"{munged_prefix}_h2"

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

        run_cmd(cmd_h2, "LDSC_HERITABILITY")

        log("\n🎉 Pipeline finished successfully!")
        log(f"Munged file: {munged_file}")
        log(f"h2 output: {out_h2}")

    except Exception as e:
        log("\n❌ Pipeline FAILED")
        log(str(e))
        raise
