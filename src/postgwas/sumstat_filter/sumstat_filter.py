import os
import subprocess
from typing import Optional
import os, subprocess, io
from pathlib import Path
from postgwas.utils.main import run_cmd




def filter_gwas_vcf_bcftools(
    vcf_path: str,
    output_folder: str,
    output_prefix: str,
    pval_cutoff: Optional[float] = None,
    maf_cutoff: Optional[float] = None,
    allelefreq_diff_cutoff: Optional[float] = None,
    info_cutoff: Optional[float] = None,
    external_af_name: str = "EUR",
    include_indels: bool = True,
    include_palindromic: bool = True,
    palindromic_af_lower: float = 0.4,
    palindromic_af_upper: float = 0.6,
    remove_mhc: bool = False,
    mhc_chrom: str = "6",
    mhc_start: int = 25000000,
    mhc_end: int = 34000000,
    threads: int = 5,
    max_mem="5G",
) -> str:
    """
    Filter VCF using bcftools — NO stdout printing.
    All logs are written to:
        output_folder/logs/<output_prefix>_filter_gwas_vcf_bcftools.log
    """
    # ============================================================
    # Create log directory + log buffer
    # ============================================================
    log_dir = Path(output_folder) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{output_prefix}_filter_gwas_vcf_bcftools.log"
    log_buffer = io.StringIO()
    def log_print(*args):
        log_buffer.write(" ".join(str(a) for a in args) + "\n")
    # ============================================================
    # Header: Start
    if not os.path.exists(vcf_path):
        msg = f"❌ ERROR: Input VCF not found: {vcf_path}"
        log_print(msg)
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())
        raise FileNotFoundError(msg)
    os.makedirs(output_folder, exist_ok=True)
    output_vcf = os.path.join(output_folder, f"{output_prefix}_filtered.vcf.gz")
    # ============================================================
    # Count variants before filtering
    # ============================================================
    pre_cmd = f"bcftools view --threads {threads} -H {vcf_path} | wc -l"
    pre = run_cmd(pre_cmd)
    try:
        pre_variants = int(pre.stdout.strip())
    except:
        pre_variants = None
    log_print(f"📊 Variants BEFORE filtering: {pre_variants}")
    log_print("")
    # ============================================================
    # Build inclusion expression
    # ============================================================
    include_parts = []
    if pval_cutoff is not None:
        include_parts.append(f"(FORMAT/LP >= {pval_cutoff})")
    if maf_cutoff is not None:
        include_parts.append(f"(FORMAT/AF >= {maf_cutoff} & FORMAT/AF <= {1 - maf_cutoff})")
    if info_cutoff is not None:
        include_parts.append(f"(FORMAT/SI >= {info_cutoff})")
    if allelefreq_diff_cutoff is not None:
        include_parts.append(
            f"(abs(FORMAT/AF - INFO/{external_af_name}) <= {allelefreq_diff_cutoff})"
        )
    include_expr = " & ".join(include_parts) if include_parts else "1"
    log_print("🔧 Include expression (bcftools -i):")
    log_print(include_expr)
    log_print("")
    # ============================================================
    # Build exclusion expression
    # ============================================================
    palindromic_logic = (
        '((REF=="A" & ALT=="T") | (REF=="T" & ALT=="A") | '
        '(REF=="C" & ALT=="G") | (REF=="G" & ALT=="C"))'
    )
    exclude_expr = None
    if not include_palindromic:
        exclude_expr = (
            f"({palindromic_logic} & "
            f"(FORMAT/AF >= {palindromic_af_lower} & FORMAT/AF <= {palindromic_af_upper}))"
        )
    log_print("🔧 Exclude expression (bcftools -e):")
    log_print(str(exclude_expr))
    log_print("")
    # ============================================================
    # Build pipeline
    # ============================================================
    cmd_parts = []
    cmd_parts.append(f"bcftools view --threads {threads} -i '{include_expr}' {vcf_path}")
    if not include_indels:
        cmd_parts.append(f"bcftools view --threads {threads} --types snps")
    if exclude_expr:
        cmd_parts.append(f"bcftools view --threads {threads} -e '{exclude_expr}'")
    if remove_mhc:
        mhc_bed = os.path.join(output_folder, f"{output_prefix}_mhc_exclude.bed")
        with open(mhc_bed, "w") as f:
            f.write(f"{mhc_chrom}\t{mhc_start}\t{mhc_end}\n")
        cmd_parts.append(f"bcftools view --threads {threads} -T ^{mhc_bed}")
    cmd_parts.append(f"bcftools sort --temp-dir {output_folder} --max-mem {max_mem}")
    cmd_parts.append("bgzip -c")
    pipeline = " | ".join(cmd_parts) + f" > {output_vcf} && tabix -f -p vcf {output_vcf}"
    log_print("🚀 Full bcftools pipeline:")
    log_print(pipeline)
    log_print("")
    # ============================================================
    # Execute pipeline
    # ============================================================
    result = run_cmd(pipeline)
    if result.returncode != 0:
        log_print("❌ Error during bcftools pipeline execution.")
        log_print("STDERR:")
        log_print(result.stderr)
        with open(log_file, "w") as f:
            f.write(log_buffer.getvalue())
        raise RuntimeError("bcftools pipeline failed")
    # ============================================================
    # Count variants after filtering
    # ============================================================
    post_cmd = f"bcftools view --threads {threads} -H {output_vcf} | wc -l"
    post = run_cmd(post_cmd)
    try:
        post_variants = int(post.stdout.strip())
    except:
        post_variants = None
    log_print(f"✅ Variants AFTER filtering: {post_variants}")
    log_print(f"💾 Filtered VCF saved to: {output_vcf}")
    log_print("\n🎉 bcftools filtering completed.")
    # ============================================================
    # Write log file
    # ============================================================
    with open(log_file, "w") as f:
        f.write(log_buffer.getvalue())
    return output_vcf

