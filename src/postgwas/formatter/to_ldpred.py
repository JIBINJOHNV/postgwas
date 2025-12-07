import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def run_ldpred_chromosome(args):
    chrom, sumstat_vcf, sample_name, output_folder = args
    chr_output = output_folder / f"{sample_name}_chr{chrom}_predld_input.tsv"
    cmd = f"""
    (
      echo "snp\\tchr\\tpos\\tA1\\tA2\\tbeta\\tSE\\tNC\\tSS\\tAF\\tLP\\tSI" && \
      bcftools view \
          -r {chrom} \
          --min-alleles 2 --max-alleles 2 "{sumstat_vcf}" | \
      bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%NC]\\t[%SS]\\t[%AF]\\t[%LP]\\t[%SI]\\n' | \
      sed 's|:|_|g'
    ) > "{chr_output}"
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return f"❌ Chr {chrom} failed:\n{result.stderr}"
    return f"✅ Chr {chrom} done: {chr_output}"



def vcf_to_ldpred(
    sumstat_vcf: str,
    output_folder: str,
    sample_name: str,
    chrom_list=range(1, 23),
    nthreads: int = 8,
):
    output_folder = Path(output_folder) / "ldpred"
    output_folder.mkdir(parents=True, exist_ok=True)
    # Prepare task arguments
    task_args = [
        (chrom, sumstat_vcf, sample_name, output_folder)
        for chrom in chrom_list
    ]
    messages = []
    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        futures = {
            executor.submit(run_ldpred_chromosome, args): args[0]
            for args in task_args
        }
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                msg = future.result()
                #print(msg)
                messages.append(msg)
            except Exception as e:
                err = f"💥 Error in chromosome {chrom}: {e}"
                print(err)
                messages.append(err)