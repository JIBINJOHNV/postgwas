import subprocess
import polars as pl
from pathlib import Path

def vcf_to_magma(sumstat_vcf: str, output_folder: str, sample_name: str):
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
        
    Returns
    -------
    dict
        Paths to the generated files.
    """
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define log file path based on sample name and function context
    log_file_path = output_dir / f"{sample_name}_magma_input.log"

    # Helper function to write to log file instead of print
    def log_message(msg: str):
        with open(log_file_path, "a", encoding="utf-8") as f:
            f.write(f"{msg}\n")

    # Clear previous log or start new session
    log_message(f"--- Starting conversion for {vcf_path.name} ---")

    # Intermediate file
    raw_output_file = output_dir / f"{sample_name}_magma_input.tsv"

    # Build bcftools pipeline
    # Note: We output raw columns and handle formatting in Polars for safety.
    # We use -e 'set -o pipefail' to catch bcftools errors even if sed succeeds.
    cmd = f"""
    set -o pipefail
    (
        printf "SNP\\tCHR\\tBP\\tREF\\tALT\\tBETA\\tSE\\tN_COL\\tAF\\tLP\\n"
        bcftools view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
        bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\t[%LP]\\n' | \
        sed 's|:|_|g'
    ) > "{raw_output_file}"
    """
    
    log_message(f"Processing {vcf_path.name}...")
    
    try:
        # executable='/bin/bash' is required for pipefail.
        # stderr is redirected to the log file to capture bcftools errors without printing to screen.
        with open(log_file_path, "a") as log_f:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stderr=log_f)
    except subprocess.CalledProcessError as e:
        log_message(f"Conversion failed for {vcf_path}: {e}")
        raise e

    # ---------- Read file ----------
    # We explicitly define null_values to handle bcftools '.' output
    # We enforce dtypes to prevent mixed-type errors on chromosomes (e.g. 1 vs X)
    try:
        df = pl.read_csv(
            raw_output_file, 
            separator="\t",
            null_values=["."],
            schema_overrides={
                "CHR": pl.String,
                "BP": pl.Int64, 
                "LP": pl.Float64,
                "N_COL": pl.Float64
            }
        )
    except Exception as e:
        error_msg = f"Failed to read intermediate file. Check VCF content. Error: {e}"
        log_message(error_msg)
        raise ValueError(error_msg)

    # ---------- Check for required columns ----------
    required_cols = ["SNP", "CHR", "BP", "REF", "ALT", "LP", "BETA", "SE", "AF", "N_COL"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        error_msg = f"Missing columns in {raw_output_file}: {', '.join(missing)}"
        log_message(error_msg)
        raise ValueError(error_msg)

    # ---------- Convert -log10(P) → P ----------
    # Protection against underflow: 
    # If LP > 320, 10^-LP becomes 0.0 in standard float64. 
    # MAGMA may reject P=0. We clamp it to 1e-320 (near min float).
    df = df.with_columns(
        pl.when(pl.col("LP") > 320)
        .then(1e-320)
        .otherwise(10 ** (-pl.col("LP")))
        .alias("P")
    )

    # ---------- Select and reorder ----------
    df_clean = df.select([
        pl.col("SNP"),
        pl.col("CHR"),
        pl.col("BP"),
        pl.col("REF"),
        pl.col("ALT"),
        pl.col("P"),
        pl.col("BETA"),
        pl.col("SE"),
        pl.col("AF"),
        pl.col("N_COL")
    ])

    # ---------- Define output files ----------
    snp_loc_file = output_dir / f"{sample_name}_magma_snp_loc.tsv"
    pval_file = output_dir / f"{sample_name}_magma_P_val.tsv"

    # ---------- Write files ----------
    # MAGMA .loc file: SNP CHR BP
    df_clean.select(["SNP", "CHR", "BP"]).write_csv(snp_loc_file, separator="\t")
    
    # MAGMA .pval file: SNP P N (N is optional but recommended if sample size varies)
    df_clean.select(["SNP", "P", "N_COL"]).write_csv(pval_file, separator="\t")

    log_message(f"MAGMA input files successfully generated in: {output_dir}")
    log_message(f"SNP location file: {snp_loc_file.name}")
    log_message(f"P-value file: {pval_file.name}")

    return {
        "snp_loc_file": str(snp_loc_file),
        "pval_file": str(pval_file),
    }