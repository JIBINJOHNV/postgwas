import subprocess
from pathlib import Path
import polars as pl




## reason to use effective sample size and sample prevalence 0.5 ; https://groups.google.com/g/ldsc_users/c/yJT-_qSh_44/m/MmKKJYsBAwAJ

def vcf_to_ldsc(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert a summary-statistics VCF into LDSC-ready tab-delimited format.
    """
    from pathlib import Path
    import subprocess

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{sample_name}_ldsc_input.tsv"
    cmd = f"""
    {{
        printf "SNP\\tA2\\tA1\\tZ\\tP\\tN\\tFRQ\\tINFO\\n"
        bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \\
        bcftools query -f '%ID\\t%REF\\t%ALT\\t[%EZ]\\t[%LP]\\t[%NEF]\\t[%AF]\\t[%SI]\\n' | \\
        sed 's|:|_|g' | \\
        awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
            # Convert -log10(P) to P
            raw_p = exp(-$5 * log(10))
            $5 = sprintf("%.6g", raw_p)
            print
        }}'
    }} > "{output_file}"
    """
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    except subprocess.CalledProcessError as e:
        print(f"❌ Conversion failed for {vcf_path} : {e}")
