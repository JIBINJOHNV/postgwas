import subprocess
from pathlib import Path
import polars as pl
import numpy as np



def ld_clump_by_regions(
    sumstat_vcf: str,
    output_folder: str,
    sample_name: str,
    population: str = "EUR",
    bcftools: str ="bcftools",
    nthreads=4
):
    """
    Extract summary statistics from a VCF and perform LD-block-based pruning
    using the population-specific LDblock INFO field (EUR_LDblock, AFR_LDblock, or EAS_LDblock),
    keeping the SNP with the highest LP (-log10 p-value) per block.

    Adds:
      - raw P-value column (10^-LP)
      - START and END columns parsed from LDblock string (e.g., "EUR-6_26791233_28017819")

    The intermediate file is compressed (.tsv.gz).

    Output
    ------
    <output_folder>/<sample_name>_vcf.tsv.gz
    <output_folder>/<sample_name>_LDpruned_<population>.tsv
    <output_folder>/<sample_name>_LDpruned_<population>_sig.tsv
    """
    population = population.upper()
    if population not in {"EUR", "AFR", "EAS"}:
        raise ValueError("Population must be one of: 'EUR', 'AFR', or 'EAS'")
    
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_out = output_dir / f"{sample_name}_vcf.tsv.gz"
    pruned_out = output_dir / f"{sample_name}_LDpruned_{population}.tsv"
    pruned_sig_out = output_dir / f"{sample_name}_LDpruned_{population}_sig.tsv"
    # --------------------- 1️⃣ Extract INFO fields ---------------------
    #print("🔹 Extracting fields from VCF with bcftools...")
    cmd = f"""
    {{
        printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tBETA\\tSE\\tNEF\\tAF\\tAFR\\tEAS\\tEUR\\tSAS\\tLP\\tNC\\tNCO\\tEUR_LDblock\\tAFR_LDblock\\tEAS_LDblock\\tBCSQ\\n"
        {bcftools} view --threads {nthreads} --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
        {bcftools} query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\t[%AFR]\\t[%EAS]\\t[%EUR]\\t[%SAS]\\t[%LP]\\t[%NC]\\t[%NCO]\\t%INFO/EUR_LDblock\\t%INFO/AFR_LDblock\\t%INFO/EAS_LDblock\\t[%BCSQ]\\n' | \
        sed 's|:|_|g'
    }} | gzip -c > "{raw_out}"
    """
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        #print(f"✅ Extracted and compressed summary stats → {raw_out}")
    except subprocess.CalledProcessError as e:
        print(f"❌ bcftools extraction failed: {e}")
        return None
    # --------------------- 2️⃣ Load with Polars ---------------------
    #print("🔹 Loading compressed TSV with Polars...")
    df = pl.read_csv(raw_out, separator="\t", null_values=[".", "NA"], infer_schema_length=20000)
    # --------------------- 3️⃣ Prepare LP, P-value, and LDblock parsing ---------------------
    ld_field = f"{population}_LDblock"
    if ld_field not in df.columns:
        raise ValueError(f"Expected LD block column '{ld_field}' not found in file.")
    df = (
        df.with_columns([
            pl.col(ld_field).alias("LDblock"),
            pl.col("LP").cast(pl.Float64),
        ])
        .filter(pl.col("LDblock").is_not_null())
        .with_columns([
            (10 ** (-pl.col("LP"))).alias("P_value"),
            # Split LDblock to extract start and end coordinates
            pl.col("LDblock").str.split("_").alias("LD_parts")
        ])
        .with_columns([
            pl.col("LD_parts").list.get(1).cast(pl.Int64).alias("START"),
            pl.col("LD_parts").list.get(2).cast(pl.Int64).alias("END"),
        ])
        .drop("LD_parts")  # keep LDblock, drop temp list
    )
    # --------------------- 4️⃣ LD-block-based pruning ---------------------
    pruned = (
        df.sort("LP", descending=True)
          .group_by("LDblock", maintain_order=True)
          .head(1) )
    # --------------------- 5️⃣ Genome-wide significant subset ---------------------
    pruned_sig = pruned.filter(pl.col("LP") > 7.3)  # p < 5×10⁻⁸
    # --------------------- 6️⃣ Save outputs ---------------------
    pruned.write_csv(pruned_out, separator="\t")
    pruned_sig.write_csv(pruned_sig_out, separator="\t")
    # print(f"✅ LD-block pruned (by LP, {population}) SNPs → {pruned_out}")
    # print(f"✅ Significant (LP > 7.3, {population}) SNPs → {pruned_sig_out}")
    # print(f"✅ Added columns: START, END (from LDblock string)")
    return str(pruned_out)




def ld_clump_standard():
    pass 



# ld_prune(
#     sumstat_vcf="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/PGC3_SCZ_european_GRCh37_merged.vcf.gz",
#     output_folder="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/PGC3_SCZ_european/",
#     sample_name="PGC3_SCZ_european",
#     population="EUR"
# )
