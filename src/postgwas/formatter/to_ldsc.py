import subprocess
from pathlib import Path
import polars as pl




## reason to use effective sample size and sample prevalence 0.5 ; https://groups.google.com/g/ldsc_users/c/yJT-_qSh_44/m/MmKKJYsBAwAJ
from pathlib import Path
import subprocess
import textwrap


def vcf_to_ldsc(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert a summary-statistics VCF into LDSC-ready tab-delimited format.

    Output
    ------
    <output_folder>/<sample_name>_ldsc_input.tsv
    <output_folder>/<sample_name>_vcf_to_ldsc.log
    """

    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{sample_name}_ldsc_input.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_ldsc.log"

    cmd = textwrap.dedent(f"""
        {{
            printf "SNP\\tA2\\tA1\\tZ\\tP\\tN\\tFRQ\\tINFO\\n"
            bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            bcftools query -f '%ID\\t%REF\\t%ALT\\t[%EZ]\\t[%LP]\\t[%NEF]\\t[%AF]\\t[%SI]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                # Convert -log10(P) to P
                raw_p = exp(-$5 * log(10))
                $5 = sprintf("%.6g", raw_p)
                print
            }}'
        }} > "{output_file}"
    """)

    try:
        with open(log_file, "w") as lf:
            lf.write("### vcf_to_ldsc command\n\n")
            lf.write(cmd + "\n")
            lf.write("\n### Command output (stdout / stderr)\n\n")

            subprocess.run(
                cmd,
                shell=True,
                executable="/bin/bash",
                check=True,
                stdout=lf,
                stderr=lf,
            )

    except subprocess.CalledProcessError as e:
        with open(log_file, "a") as lf:
            lf.write("\n❌ Conversion failed\n")
            lf.write(str(e) + "\n")

        raise RuntimeError(
            f"vcf_to_ldsc failed for {vcf_path}. "
            f"See log file: {log_file}"
        )

    return {
        "ldsc_file": str(output_file),
        "log_file": str(log_file),
    }
