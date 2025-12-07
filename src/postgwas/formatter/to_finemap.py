import subprocess
from pathlib import Path
import polars as pl



def vcf_to_finemap(sumstat_vcf: str, output_folder: str, sample_name: str):
    """
    Convert a single summary-statistics VCF file into a MAGMA-compatible
    tab-delimited text file.
    Parameters
    ----------
    sumstat_vcf : str
        Path to the input VCF (e.g., "PGC3_SCZ_european_chr1.vcf.gz")
    output_folder : str
        Folder where the converted text file will be written
    sample_name : str
        Prefix for output naming (e.g., "PGC3_SCZ_european")
    Output
    ------
    <output_folder>/<sample_name>_<chromosome>_magma_input.txt
    """
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Output file naming
    output_file = output_dir / f"{sample_name}_finemap.tsv"
    # Build bcftools → query → sed pipeline
    cmd = f"""
    {{
        printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tEZ\\tLP\\tNEF\\tAF\\n"
        bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
        bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%EZ]\\t[%LP]\\t[%NEF]\\t[%AF]\\n' | \
        sed 's|:|_|g'
    }} > "{output_file}"
    """
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    except subprocess.CalledProcessError as e:
        print(f"❌ Conversion failed for {vcf_path} : {e}")


