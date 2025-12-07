
import os,subprocess,re,sys,io
from typing import Optional
import textwrap
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


def run_bcftools_munge(
    output_dir: str,
    gwas_outputname: str,
    chr: str,
    genome_fasta_file: str,
    grch_version:str,
) -> str:
    """
    Run bcftools +munge → bcftools norm → bcftools sort to convert GWAS summary stats to VCF.

    Parameters
    ----------
    output_dir : str
        Directory where input/output files are located.
    gwas_outputname : str
        Base GWAS dataset name.
    chr : str
        Chromosome number or label.
    genome_fasta_file : str
        Path to genome FASTA reference file.

    Returns
    -------
    str
        Path to the final VCF.gz file.
    """
    # Paths
    input_file = f"{output_dir}/{gwas_outputname}_chr{chr}_vcf_input.tsv"
    column_file = f"{output_dir}/{gwas_outputname}_chr{chr}_column_mapping.tsv"
    output_vcf = f"{output_dir}/{gwas_outputname}_chr{chr}_{grch_version}.vcf.gz"
    
    # Ensure directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the bcftools command
    cmd = (
        f"bcftools +munge --no-version "
        f"--columns-file {column_file} "
        f"--fasta-ref {genome_fasta_file} "
        f"--sample-name {gwas_outputname} "
        f"{input_file} | "
        f"bcftools norm -m-any -d exact | "
        f"bcftools sort -Oz -o {output_vcf} --write-index=tbi "
    )
    print(f"🚀 Running BCFtools command:\n{cmd}\n")
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print(f"✅ Completed: {output_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running bcftools pipeline:\n{e}")
        return "Error"
    return output_vcf



def run_bcftools_annot(
    output_dir: str,
    gwas_outputname: str,
    chromosome: str,
    external_eaf_file: str,
    default_dbsnp_file: str,
    genome_fasta_file: str,
    target_genome_fasta_file: str,
    gff_file: str,
    chain_file: str,
    grch_version: str,
    threads: int = 5
) -> str:
    # -------------------------------------------------------
    # 0. CREATE LOG FILE
    # -------------------------------------------------------
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_chr{chromosome}_bcftools_annot.log"
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    # -------------------------------------------------------
    # 1. ALL LOGGING INSIDE THE with BLOCK
    # -------------------------------------------------------
    try:
        with open(log_file, "w") as f:
            # Redirect process output into the open file
            sys.stdout = f
            sys.stderr = f
            try:
                print(f"▶ Starting bcftools annotation for chromosome {chromosome}")
                # -------------------------------------------------------
                # Paths
                # -------------------------------------------------------
                outdir = output_dir
                os.makedirs(outdir, exist_ok=True)
                input_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{grch_version}.vcf.gz"
                output1_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{grch_version}_ID.vcf.gz"
                output2_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{grch_version}_EAF.vcf.gz"
                output3_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{grch_version}_CSQ.vcf.gz"
                # -------------------------------------------------------
                # Step 1: dbSNP annotation
                # -------------------------------------------------------
                print("🧬 Step 1/4: dbSNP annotation...")
                subprocess.run(
                    [
                        "bcftools", "annotate",
                        "--threads", str(threads),
                        "--annotations", default_dbsnp_file,
                        "--columns", "CHROM,POS,REF,ALT,ID",
                        "--output-type", "z",
                        "--output", output1_vcf,
                        "--write-index=tbi",
                        input_vcf
                    ],
                    stdout=f, stderr=f, check=True
                )
                # -------------------------------------------------------
                # Step 2: AF annotation
                # -------------------------------------------------------
                print("🌍 Step 2/4: AF annotation...")
                subprocess.run(
                    [
                        "bcftools", "annotate",
                        "--threads", str(threads),
                        "--annotations", external_eaf_file,
                        "--columns", "CHROM,POS,REF,ALT,INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS",
                        "--output-type", "z",
                        "--output", output2_vcf,
                        "--write-index=tbi",
                        output1_vcf
                    ],
                    stdout=f, stderr=f, check=True
                )
                # -------------------------------------------------------
                # Step 3: bcftools csq
                # -------------------------------------------------------
                print("🧫 Step 3/4: Functional consequence annotation (bcftools csq)...")
                subprocess.run(
                    [
                        "bcftools", "csq",
                        "--fasta-ref", genome_fasta_file,
                        "--gff-annot", gff_file,
                        "--threads", str(threads),
                        "--unify-chr-names", "-,chr,-",
                        "--output-type", "z",
                        "--output", input_vcf,
                        "--write-index=tbi",
                        output2_vcf
                    ],
                    stdout=f, stderr=f, check=True
                )
                # -------------------------------------------------------
                # Determine target genome build
                # -------------------------------------------------------
                if grch_version == "GRCh37":
                    target_build = "GRCh38"
                else:
                    target_build = "GRCh37"
                target_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{target_build}.vcf.gz"
                reject_vcf = f"{outdir}/{gwas_outputname}_chr{chromosome}_{target_build}_notlifted.vcf.gz"
                # -------------------------------------------------------
                # Step 4: Liftover
                # -------------------------------------------------------
                print(f"🧭 Step 4/4: liftover {grch_version} → {target_build}...")
                liftover_cmd = f"""
                    bcftools +liftover --no-version \
                        --output-type u {input_vcf} -- \
                        --src-fasta-ref {genome_fasta_file} \
                        --fasta-ref {target_genome_fasta_file} \
                        --chain {chain_file} \
                        --reject {reject_vcf} \
                        --reject-type z \
                    | bcftools view -e 'INFO/SWAP ==1 || INFO/SWAP ==-1' | \
                        bcftools sort \
                        --output {target_vcf} \
                        --output-type z \
                        --write-index=tbi
                """
                subprocess.run(["bash", "-c", liftover_cmd], stdout=f, stderr=f, check=True)
                print(f"✅ DONE → {target_vcf}")
            except Exception as pipeline_error:
                f.write(f"\n❌ ERROR inside run_bcftools_annot (chr{chromosome}): {pipeline_error}\n")
    finally:
        # Always restore stdout/stderr after the file is closed.
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        f.close() 
    return f"Completed bcftools annotation for chr{chromosome}"



import io
import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def concat_vcfs_by_build(output_dir: str, gwas_outputname: str, mode: str = "concurrent"):
    """
    Concatenate chromosome-wise VCFs for GRCh37/38 and notlifted files.
    All logs go into logs/{sample}_concat_vcfs_by_build.log.
    Works on Linux + macOS.
    """

    # ------------------------------------------------------------------
    # Ensure Path object
    # ------------------------------------------------------------------
    outdir = Path(output_dir).resolve()
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{gwas_outputname}_concat_vcfs_by_build.log"

    # Log buffer (avoid verbose screen printing)
    log_buffer = io.StringIO()

    def log_write(*args):
        msg = " ".join(str(a) for a in args)
        log_buffer.write(msg + "\n")

    # ------------------------------------------------------------------
    # Helper: chromosome sorting key
    # ------------------------------------------------------------------
    def chr_sort_key(path: Path):
        m = re.search(r"chr(\d+|X|Y|MT)", path.stem)
        if not m:
            return 999
        v = m.group(1)
        if v == "X": return 23
        if v == "Y": return 24
        if v == "MT": return 25
        return int(v)

    # ------------------------------------------------------------------
    # Collect chromosome VCFs
    # ------------------------------------------------------------------
    all_vcfs = list(outdir.glob(f"{gwas_outputname}_chr*_*GRCh*.vcf.gz"))
    if not all_vcfs:
        log_write("❌ No chromosome VCF files found in output directory.")
        with open(log_file, "w") as f: f.write(log_buffer.getvalue())
        raise FileNotFoundError("No chromosome VCFs found")

    valid_vcfs = [
        f for f in all_vcfs
        if "_notlifted" not in f.name and
           (f.name.endswith("_GRCh37.vcf.gz") or f.name.endswith("_GRCh38.vcf.gz"))
    ]
    notlifted_vcfs = [f for f in all_vcfs if "_notlifted" in f.name]

    # Split by genome build
    grch37_vcfs = sorted([p for p in valid_vcfs if "GRCh37" in p.name], key=chr_sort_key)
    grch38_vcfs = sorted([p for p in valid_vcfs if "GRCh38" in p.name], key=chr_sort_key)
    notlifted_37 = sorted([p for p in notlifted_vcfs if "GRCh37" in p.name], key=chr_sort_key)
    notlifted_38 = sorted([p for p in notlifted_vcfs if "GRCh38" in p.name], key=chr_sort_key)

    merged_paths = {
        "grch37":     outdir / f"{gwas_outputname}_GRCh37_merged.vcf.gz",
        "grch38":     outdir / f"{gwas_outputname}_GRCh38_merged.vcf.gz",
        "notlifted_37": outdir / f"{gwas_outputname}_GRCh37_notlifted_merged.vcf.gz",
        "notlifted_38": outdir / f"{gwas_outputname}_GRCh38_notlifted_merged.vcf.gz",
    }

    # ------------------------------------------------------------------
    # bcftools concat runner
    # ------------------------------------------------------------------
    def merge_vcfs(tag: str, vcf_list, out_path: Path):
        if not vcf_list:
            log_write(f"⚠️ No VCF files found for {tag}. Skipping.")
            return None

        log_write(f"🔧 [{tag}] Running bcftools concat on {len(vcf_list)} files…")

        result = subprocess.run(
            [
                "bcftools", "concat",
                "--output-type", "z",
                "--output", str(out_path),
                "--threads", "4",
                "--write-index=tbi",
                *[str(f) for f in vcf_list],
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Log bcftools messages
        if result.stdout.strip():
            log_write(f"[bcftools stdout:{tag}] {result.stdout.strip()}")
        if result.stderr.strip():
            log_write(f"[bcftools stderr:{tag}] {result.stderr.strip()}")

        if result.returncode != 0:
            log_write(f"❌ bcftools concat failed for {tag} (exit {result.returncode})")
            return None

        log_write(f"✅ [{tag}] Merge completed → {out_path}")

        # Clean chromosome-level files
        for f in vcf_list:
            try:
                f.unlink(missing_ok=True)
                tbi = f.with_suffix(f.suffix + ".tbi")
                if tbi.exists():
                    tbi.unlink()
            except Exception:
                pass

        log_write(f"🧹 [{tag}] Cleaned per-chromosome VCF + index files.")
        return str(out_path)

    # ------------------------------------------------------------------
    # Execute merges
    # ------------------------------------------------------------------
    tasks = []

    if mode == "concurrent":
        with ThreadPoolExecutor(max_workers=4) as ex:
            if grch37_vcfs:
                tasks.append(ex.submit(merge_vcfs, "GRCh37", grch37_vcfs, merged_paths["grch37"]))
            if grch38_vcfs:
                tasks.append(ex.submit(merge_vcfs, "GRCh38", grch38_vcfs, merged_paths["grch38"]))
            if notlifted_37:
                tasks.append(ex.submit(merge_vcfs, "GRCh37_notlifted", notlifted_37, merged_paths["notlifted_37"]))
            if notlifted_38:
                tasks.append(ex.submit(merge_vcfs, "GRCh38_notlifted", notlifted_38, merged_paths["notlifted_38"]))

            results = [t.result() for t in tasks]
    else:
        results = []
        if grch37_vcfs:
            results.append(merge_vcfs("GRCh37", grch37_vcfs, merged_paths["grch37"]))
        if grch38_vcfs:
            results.append(merge_vcfs("GRCh38", grch38_vcfs, merged_paths["grch38"]))
        if notlifted_37:
            results.append(merge_vcfs("GRCh37_notlifted", notlifted_37, merged_paths["notlifted_37"]))
        if notlifted_38:
            results.append(merge_vcfs("GRCh38_notlifted", notlifted_38, merged_paths["notlifted_38"]))

    # ------------------------------------------------------------------
    # Cleanup: remove any leftover chr*.vcf.gz or .tbi
    # ------------------------------------------------------------------
    for f in outdir.glob(f"{gwas_outputname}_chr*.vcf.gz*"):
        try:
            f.unlink()
        except Exception:
            pass

    log_write("🎯 All merges completed successfully.")

    # Write logfile
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())

    return merged_paths








# def concat_vcfs_by_build(output_dir: str, gwas_outputname: str, mode: str = "concurrent"):
#     """
#     Concatenate chromosome-wise VCFs and delete intermediate files.
#     All messages are written to a logfile only (no screen printing).
#     """

#     # ------------------------------------------------------------------
#     # Setup log file
#     # ------------------------------------------------------------------
#     outdir = output_dir
#     log_dir = Path(output_dir) / "logs"
#     log_dir.mkdir(parents=True, exist_ok=True)
#     log_file = log_dir / f"{gwas_outputname}_concat_vcfs_by_build.log"

#     log_buffer = io.StringIO()

#     def log_write(*args):
#         msg = " ".join(str(a) for a in args)
#         log_buffer.write(msg + "\n")

#     # ------------------------------------------------------------------
#     # Helper sorting function
#     # ------------------------------------------------------------------
#     def chr_sort_key(path):
#         m = re.search(r"chr(\d+|X|Y|MT)", path.stem)
#         if not m:
#             return 999
#         val = m.group(1)
#         if val == "X": return 23
#         if val == "Y": return 24
#         if val == "MT": return 25
#         return int(val)

#     # ------------------------------------------------------------------
#     # Collect chromosome-level VCFs
#     # ------------------------------------------------------------------
#     all_vcfs = list(outdir.glob(f"{gwas_outputname}_chr*_*GRCh*.vcf.gz"))
#     if not all_vcfs:
#         log_write("❌ No chromosome VCF files found in output directory.")
#         with open(log_file, "w") as f:
#             f.write(log_buffer.getvalue())
#         raise FileNotFoundError("No chromosome VCFs found")

#     valid_vcfs = [
#         f for f in all_vcfs
#         if "_notlifted" not in f.name
#         and any(f.name.endswith(suffix) for suffix in ["_GRCh37.vcf.gz", "_GRCh38.vcf.gz"])
#     ]

#     notlifted_vcfs = [f for f in all_vcfs if "_notlifted" in f.name]

#     # Separate by build
#     grch37_vcfs = sorted([p for p in valid_vcfs if "GRCh37" in p.name], key=chr_sort_key)
#     grch38_vcfs = sorted([p for p in valid_vcfs if "GRCh38" in p.name], key=chr_sort_key)
#     notlifted_37 = sorted([p for p in notlifted_vcfs if "GRCh37" in p.name], key=chr_sort_key)
#     notlifted_38 = sorted([p for p in notlifted_vcfs if "GRCh38" in p.name], key=chr_sort_key)

#     merged_files = {
#         "grch37": outdir / f"{gwas_outputname}_GRCh37_merged.vcf.gz",
#         "grch38": outdir / f"{gwas_outputname}_GRCh38_merged.vcf.gz",
#         "notlifted_37": outdir / f"{gwas_outputname}_GRCh37_notlifted_merged.vcf.gz",
#         "notlifted_38": outdir / f"{gwas_outputname}_GRCh38_notlifted_merged.vcf.gz",
#     }

#     # ------------------------------------------------------------------
#     # Merge helper
#     # ------------------------------------------------------------------
#     def merge_vcfs(tag, vcf_list, out_path):
#         if not vcf_list:
#             log_write(f"⚠️ No files found for {tag}, skipping.")
#             return None

#         # Run bcftools fully silently
#         result = subprocess.run(
#             [
#                 "bcftools", "concat",
#                 "--output-type", "z",
#                 "--output", str(out_path),
#                 "--threads", "4",
#                 "--write-index=tbi",
#                 *map(str, vcf_list)
#             ],
#             stdout=subprocess.PIPE,     # silence normal output
#             stderr=subprocess.PIPE,     # silence error output
#             text=True                   # decode bytes → str
#         )
#         log_write(f"✅ [{tag}] Merged {len(vcf_list)} chromosomes → {out_path}")
#         # Log bcftools output (stdout + stderr)
#         if result.stdout:
#             log_write(f"[bcftools stdout:{tag}] {result.stdout.strip()}")
#         if result.stderr:
#             log_write(f"[bcftools stderr:{tag}] {result.stderr.strip()}")

#         if result.returncode != 0:
#             log_write(f"❌ bcftools concat failed for {tag} (exit {result.returncode})")
#             return None
    
#         # Clean up per-chromosome VCFs + tbi
#         for f in vcf_list:
#             try:
#                 f.unlink()
#                 tbi = f.with_suffix(f.suffix + ".tbi")
#                 if tbi.exists():
#                     tbi.unlink()
#             except Exception:
#                 pass

#         log_write(f"🧹 [{tag}] Cleaned {len(vcf_list)} per-chromosome VCFs.")
#         return str(out_path)

#     # ------------------------------------------------------------------
#     # Run merges concurrently or sequentially
#     # ------------------------------------------------------------------
#     results = []

#     if mode == "concurrent":
#         with ThreadPoolExecutor(max_workers=4) as executor:
#             futures = []
#             if grch37_vcfs:
#                 futures.append(executor.submit(merge_vcfs, "GRCh37", grch37_vcfs, merged_files["grch37"]))
#             if grch38_vcfs:
#                 futures.append(executor.submit(merge_vcfs, "GRCh38", grch38_vcfs, merged_files["grch38"]))
#             if notlifted_37:
#                 futures.append(executor.submit(merge_vcfs, "GRCh37_notlifted", notlifted_37, merged_files["notlifted_37"]))
#             if notlifted_38:
#                 futures.append(executor.submit(merge_vcfs, "GRCh38_notlifted", notlifted_38, merged_files["notlifted_38"]))

#             results = [f.result() for f in futures]

#     else:
#         if grch37_vcfs:
#             results.append(merge_vcfs("GRCh37", grch37_vcfs, merged_files["grch37"]))
#         if grch38_vcfs:
#             results.append(merge_vcfs("GRCh38", grch38_vcfs, merged_files["grch38"]))
#         if notlifted_37:
#             results.append(merge_vcfs("GRCh37_notlifted", notlifted_37, merged_files["notlifted_37"]))
#         if notlifted_38:
#             results.append(merge_vcfs("GRCh38_notlifted", notlifted_38, merged_files["notlifted_38"]))

#     # ------------------------------------------------------------------
#     # Extra cleanup: any leftover chr*.vcf.gz or .tbi
#     # ------------------------------------------------------------------
#     for f in outdir.glob(f"{gwas_outputname}_chr*.vcf.gz*"):
#         try:
#             f.unlink()
#         except Exception:
#             pass

#     log_write("🎯 All concatenations completed.")

#     # Save log file
#     with open(log_file, "w") as f:
#         f.write(log_buffer.getvalue())

#     return merged_files
