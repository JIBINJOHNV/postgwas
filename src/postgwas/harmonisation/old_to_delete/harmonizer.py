
import os,sys
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
import polars as pl
import textwrap
from typing import Optional, Tuple, List, Dict

from postgwas.harmonisation.io import read_sumstats, read_config,find_resource_file_path
from postgwas.harmonisation.validator import validate_gwas_config
from postgwas.harmonisation.preprocess import split_chr_pos
from postgwas.harmonisation.chr_pos_process import fix_chr_pos_column
from postgwas.harmonisation.genome_build import genome_build
from postgwas.harmonisation.qc_utils import check_eaf_or_maf, compute_effective_n
from postgwas.harmonisation.parallel import parallel_split_and_run
from postgwas.harmonisation.reference import load_reference_build
#from postgwas.harmonization.vcf_converter import run_gwas2vcf
from postgwas.harmonisation.eaf_maf import add_or_calculate_eaf
from postgwas.harmonisation.case_control import harmonize_sample_sizes 
from postgwas.harmonisation.beta_or_or import is_beta_or_or
from postgwas.harmonisation.beta_and_se_from_z import calculate_beta_and_se_from_z
from postgwas.harmonisation.calculate_z_from_beta_se import calculate_z_from_beta_se
from postgwas.harmonisation.detect_pvalue_type_process import detect_and_convert_pval,detect_pval_type,convert_pval_to_mlogp
from postgwas.harmonisation.calculate_se_from_beta_pvalue import calculate_se_from_beta_pvalue
from postgwas.harmonisation.add_or_calculate_info import add_or_calculate_info
#from postgwas.harmonization.export_gwas_sumstat import export_gwas_sumstat
from postgwas.harmonisation.add_or_fix_snp_column import add_or_fix_snp_column
from postgwas.harmonisation.sumstat_to_vcf import run_bcftools_munge,run_bcftools_annot,concat_vcfs_by_build
from postgwas.harmonisation.export_gwas_sumstat_forgwas2vcf import export_gwas_sumstat
from postgwas.harmonisation.gwastovcf_gwas2vcf import gwastovcf
from postgwas.harmonisation.clean_intermediate import clean_intermediate_files
from postgwas.config import load_config
from postgwas.harmonisation.postgwas_qc import postgwas_qc

# ===============================================================
# 1️⃣ Resource Map Builder
# ===============================================================

def build_resource_map(
    chromosome: str,
    grch_version: str,
    resource_folder: str,
    user_eaf_file: str,
    default_eaf_file: str,
    default_comparison_af_file: str,
    user_info_file: Optional[str] = "NA",
    default_info_file: Optional[str] = "NA",
    dbsnp: str = "dbSNP",
    user_eaf_column: Optional[str] = "NA",
    default_eaf_column: Optional[str] = "NA",
    default_comparison_af_column: Optional[str] = "NA",
    user_info_column: Optional[str] = "NA",
    default_info_column: Optional[str] = "NA",
) -> Dict[str, Any]:
    """
    Construct all resource file paths for one chromosome.
    User-provided EAF/INFO files override defaults if provided.
    """
    # --- Allele frequency (EAF) resources ---
    user_eaf_path = find_resource_file_path(
        input_file=user_eaf_file, resource_folder=resource_folder,
        grch_version=grch_version,chromosome=chromosome )
    
    default_eaf_path = (
        f"{resource_folder}/{grch_version}/default_af/tab_files/"
        f"{grch_version}_{default_eaf_file}_freq_chr{chromosome}.tsv.gz" )
    
    default_comparison_af_path = (
        f"{resource_folder}/{grch_version}/external_af/vcf_files/"
        f"{grch_version}_{default_comparison_af_file}_freq_chr{chromosome}.vcf.gz" )
    
    # --- INFO score resources ---
    user_info_path = find_resource_file_path(
        input_file=user_info_file, resource_folder=resource_folder,
        grch_version=grch_version,chromosome=chromosome )
    
    default_info_path = (
        f"{resource_folder}/{grch_version}/default_infoscore/tab_files/"
        f"{grch_version}_{default_info_file}_infoscore_chr{chromosome}.tsv.gz" )
    # --- Core genome references ---
    genome_fasta_path = f"{resource_folder}/{grch_version}/fasta_files/{grch_version}_chr{chromosome}.fa"
    dbsnp_path = f"{resource_folder}/{grch_version}/dbSNP/vcf_files/{grch_version}_{dbsnp}_chr{chromosome}.vcf.gz"
    annot_path = f"{resource_folder}/{grch_version}/gff_files/{grch_version}_ensembl.gff3.gz"
    # --- Chain and target FASTA for liftover ---
    if grch_version == "GRCh37":
        target_fasta = f"{resource_folder}/GRCh38/fasta_files/GRCh38_chr{chromosome}.fa"
        chain_file = f"{resource_folder}/chain_files/GRCh37_to_GRCh38.chain"
    else:
        target_fasta = f"{resource_folder}/GRCh37/fasta_files/GRCh37_chr{chromosome}.fa"
        chain_file = f"{resource_folder}/chain_files/GRCh38_to_GRCh37.chain"
    
    return {
        # EAF
        "user_eaf_file": user_eaf_path,
        "default_eaf_file": default_eaf_path,
        "default_comparison_af_file": default_comparison_af_path,
        "user_eaf_column": user_eaf_column,
        "default_eaf_column": default_eaf_column,
        "default_comparison_af_column": default_comparison_af_column,

        # INFO
        "user_info_file": user_info_path,
        "default_info_file": default_info_path,
        "user_info_column": user_info_column,
        "default_info_column": default_info_column,

        # Genome, dbSNP, chain
        "genome_fasta_file": genome_fasta_path,
        "dbsnp_file": dbsnp_path,
        "annot_path": annot_path,
        "target_fasta": target_fasta,
        "chain_file": chain_file
    }

# ===============================================================
# 2️⃣ Per-Chromosome Harmonization Pipeline
# ===============================================================

def process_one_chromosome(
    chromosome: str,
    chr_file: str,
    sample_column_dict: Dict[str, Any],
    resource_folder: str,
    grch_version: str,
    user_eaf_file: str,
    default_eaf_file: str,
    default_comparison_af_file: str,
    user_info_file: Optional[str] = "NA",
    default_info_file: Optional[str] = "NA",
    user_eaf_column: Optional[str] = "NA",
    default_eaf_column: Optional[str] = "NA",
    default_comparison_af_column: Optional[str] = "NA",
    user_info_column: Optional[str] = "NA",
    default_info_column: Optional[str] = "NA",
    dbsnp: str = "dbSNP",
    output_dir: str = ".",
    threads: int = 5,
    gwastovcf_main_script_path: str = "/app/main.py",
) -> Tuple[str, Dict[str, Any]]:
    """
    Run full harmonisation + GWAS→VCF + bcftools annotation for a single chromosome.
    Returns:
        (chromosome, qc_dict)
    """
    # -------------------------------------------------------
    # 0. SET UP LOG FILE PER CHROMOSOME
    # -------------------------------------------------------
    sample_id = sample_column_dict["gwas_outputname"]
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{sample_id}_chr{chromosome}.log"
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    try:
        with open(log_file, "w") as f:
            # Redirect
            sys.stdout = f
            sys.stderr = f
            print(f"\n🚀 Starting harmonization pipeline for chromosome {chromosome}...\n")
            # ------------------------------
            # 1. Load per-chromosome file
            # ------------------------------
            df = pl.read_csv(chr_file, separator="\t")
            df = fix_chr_pos_column(df=df, sample_column_dict=sample_column_dict, drop_mt=True)
            # ------------------------------
            # 2. Build resource map
            # ------------------------------
            res = build_resource_map(
                chromosome=chromosome,
                grch_version=grch_version,
                resource_folder=resource_folder,
                user_eaf_file=user_eaf_file,
                default_eaf_file=default_eaf_file,
                default_comparison_af_file=default_comparison_af_file,
                user_info_file=user_info_file,
                default_info_file=default_info_file,
                user_eaf_column=user_eaf_column,
                default_eaf_column=default_eaf_column,
                default_comparison_af_column=default_comparison_af_column,
                user_info_column=user_info_column,
                default_info_column=default_info_column,
                dbsnp=dbsnp,
            )
            # ------------------------------
            # 3. Harmonization steps
            # ------------------------------
            df, eaf_qc, sample_column_dict = add_or_calculate_eaf(
                df=df,
                sample_column_dict=sample_column_dict,
                eaffile=res["user_eaf_file"],
                default_eaf_file=res["default_eaf_file"],
                default_eaf_eafcolumn=res["default_eaf_column"],
            )
            df, size_qc, sample_column_dict = harmonize_sample_sizes(df, sample_column_dict)
            df, beta_qc, sample_column_dict = calculate_beta_and_se_from_z(df, sample_column_dict)
            df, or_qc, sample_column_dict = is_beta_or_or(df, sample_column_dict)
            df, pval_qc, sample_column_dict = detect_and_convert_pval(df, sample_column_dict)
            df, se_qc, sample_column_dict = calculate_se_from_beta_pvalue(df, sample_column_dict)
            df, z_qc, sample_column_dict = calculate_z_from_beta_se(df, sample_column_dict)
            df, info_qc, sample_column_dict = add_or_calculate_info(
                df=df,
                sample_column_dict=sample_column_dict,
                default_info_file=res["default_info_file"],
                info_file=res["user_info_file"],
                default_info_column=res["default_info_column"],
                info_column=res["user_info_column"],
            )
            df, sample_column_dict = add_or_fix_snp_column(df, sample_column_dict)
            # ------------------------------
            # 4. Export cleaned GWAS sumstats
            # ------------------------------
            export_gwas_summary_df = export_gwas_sumstat(
                df=df,
                sample_column_dict=sample_column_dict,
                output_dir=output_dir,
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                genome_build=grch_version,
            )
            # ------------------------------
            # 5. Convert GWAS to VCF
            # ------------------------------
            gwastovcf_command_str, gwastovcf_exit_code = gwastovcf(
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                grch_version=grch_version,
                output_folder=output_dir,
                fasta=res["genome_fasta_file"],
                dbsnp=res["dbsnp_file"],
                aliasfile="NA",
                main_script_path=gwastovcf_main_script_path,
            )
            # ------------------------------
            # 6. Annotate VCF (bcftools)
            # ------------------------------
            run_bcftools_annot(
                output_dir=output_dir,
                gwas_outputname=sample_column_dict["gwas_outputname"],
                chromosome=chromosome,
                external_eaf_file=res["default_comparison_af_file"],
                default_dbsnp_file=res["dbsnp_file"],
                genome_fasta_file=res["genome_fasta_file"],
                target_genome_fasta_file=res["target_fasta"],
                gff_file=res["annot_path"],
                chain_file=res["chain_file"],
                grch_version=grch_version,
                threads=threads,
            )
            # ------------------------------
            # 7. QC summary
            # ------------------------------
            qc_dict: Dict[str, Any] = {
                "eaf_qc": eaf_qc,
                "sample_size_qc": size_qc,
                "beta_qc": beta_qc,
                "or_qc": or_qc,
                "pval_qc": pval_qc,
                "se_qc": se_qc,
                "z_qc": z_qc,
                "info_qc": info_qc,
                "gwastovcf_exit_code": gwastovcf_exit_code,
            }
            print(f"\n✅ Completed chromosome {chromosome}.\n")
            # ⭐ FIX — Restore stdout BEFORE leaving the with-block
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr
    except Exception as e:
        # Ensure stdout is restored before printing error
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        qc_dict = {"error": str(e)}
    finally:
        # NOTHING to restore; file was auto-closed by 'with'
        pass
    return chromosome, qc_dict



# ----------------------------------------------------------------------
# 4. MAIN PARALLEL DRIVER
# ----------------------------------------------------------------------

def gwas_to_vcf_parallel(
    sumstat_file: str,
    sample_column_dict: Dict[str, Any],
    output_dir: str,
    resource_folder: str,
    default_eaf_file: str = "1000G",
    default_info_file: str = "1000G",
    default_comparison_af_file: str = "1000G",
    dbsnp: str = "dbsnp155",
    required_columns_in_sumstat: Optional[List[str]] = None,
    optional_columns_in_sumstat: Optional[List[str]] = None,
    max_workers: int = 5,
    *,
    user_eaf_file: Optional[str] = "NA",
    user_eaf_column: Optional[str] = "EUR",
    user_info_file: Optional[str] = "NA",
    user_info_column: Optional[str] = "INFO",
    default_eaf_column: str = "EUR",
    default_comparison_af_column: str = "EUR",
    default_info_column: str = "INFO",
    gwastovcf_main_script_path: str = "/app/main.py",
    grch37_file=None,
    grch38_file=None
) -> Dict[str, Dict[str, Any]]:
    """
    Full GWAS-to-VCF pipeline:
      1. Reads the summary statistics.
      2. Infers genome build.
      3. Splits by chromosome.
      4. Runs per-chromosome harmonization + VCF export in parallel.
      5. Concatenates per-chrom VCFs and cleans intermediates.

    Returns:
        per_chr_qc: dict
            {
              "total_variant_infile": ...,
              "total_variant_read": ...,
              "1": {qc dict or error},
              "2": { ... },
              ...
            }
    """
    # ------------------------------------------------------------------
    # 0. Normalize paths and create output directory
    # ------------------------------------------------------------------
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # ------------------------------------------------------------------
    # 1. Load summary stats and validate config
    # ------------------------------------------------------------------
    df, file_cvariant_count, polars_rows = read_sumstats(
        sumstat_file=sumstat_file,
        output_dir=str(output_dir)
    )
    validate_gwas_config(sample_column_dict, df)
    # ------------------------------------------------------------------
    # 2. Infer genome build
    # ------------------------------------------------------------------
    genome_build_info = genome_build(df, grch37_file, grch38_file, sample_column_dict)
    grch_version = genome_build_info["inferred_build"]
    # ------------------------------------------------------------------
    # 3. Split per chromosome (filter to autosomes + chrX)
    # ------------------------------------------------------------------
    chr_files = split_chr_pos(df, sample_gwas_dict=sample_column_dict)
    chr_files = [
        f for f in chr_files
        if "chr" in Path(f).stem.lower()
        and not any(x in Path(f).stem.lower() for x in ["chrmt", "chry"])
    ]
    print(f"🧩 Found {len(chr_files)} chromosome files to process.\n")
    # ------------------------------------------------------------------
    # 4. Init QC container
    # ------------------------------------------------------------------
    per_chr_qc: Dict[str, Dict[str, Any]] = {}
    required_columns_in_sumstat = required_columns_in_sumstat or [
        "CHR", "BP", "A1", "A2", "SNP", "BETA", "SE", "Z", "LP", "N_CON", "NEFF"
    ]
    optional_columns_in_sumstat = optional_columns_in_sumstat or ["N_CAS", "INFO"]
    per_chr_qc["total_variant_infile"] = file_cvariant_count
    per_chr_qc["total_variant_read"] = polars_rows
    # ------------------------------------------------------------------
    # 5. Parallel per-chromosome processing
    # ------------------------------------------------------------------
    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chr: Dict[Any, str] = {}
            # Submit jobs
            for chr_file in chr_files:
                chrom = Path(chr_file).stem.split("_")[-1].replace("chr", "")
                future = executor.submit(
                    process_one_chromosome,
                    chromosome=chrom,
                    chr_file=chr_file,
                    sample_column_dict=sample_column_dict.copy(),
                    resource_folder=resource_folder,
                    grch_version=grch_version,
                    user_eaf_file=user_eaf_file,
                    default_eaf_file=default_eaf_file,
                    default_comparison_af_file=default_comparison_af_file,
                    user_info_file=user_info_file,
                    default_info_file=default_info_file,
                    user_eaf_column=user_eaf_column,
                    default_eaf_column=default_eaf_column,
                    default_comparison_af_column=default_comparison_af_column,
                    user_info_column=user_info_column,
                    default_info_column=default_info_column,
                    dbsnp=dbsnp,
                    output_dir=str(output_dir),
                    gwastovcf_main_script_path=gwastovcf_main_script_path,
                )
                future_to_chr[future] = chrom
            # Collect results
            for future in as_completed(future_to_chr):
                chrom = future_to_chr[future]
                exc = future.exception()
                if exc is not None:
                    # Worker failed – capture error but DO NOT re-raise
                    print(f"❌ [ERROR] Chromosome {chrom} failed: {exc}")
                    per_chr_qc[chrom] = {"error": str(exc)}
                    continue
                # Safe to get result now
                chrom_res, qc_dict = future.result()
                per_chr_qc[chrom_res] = qc_dict
                print(f"✅ Completed chromosome {chrom_res}")
    finally:
        # ------------------------------------------------------------------
        # 6. Always run concatenation + cleanup, even if some chromosomes fail
        # ------------------------------------------------------------------
        print("\n📦 Running concatenation and cleanup steps...\n")
        try:
            concat_vcfs_by_build(
                output_dir=output_dir,
                gwas_outputname=sample_column_dict["gwas_outputname"],
                mode="concurrent"
            )
        except Exception as e:
            print(f"⚠️ [WARN] concat_vcfs_by_build failed: {e}")
        try:
            clean_intermediate_files(
                output_dir=output_dir,
                gwas_outputname=sample_column_dict["gwas_outputname"]
            )
        except Exception as e:
            print(f"⚠️ [WARN] clean_intermediate_files failed: {e}")
        print("\n🏁 Pipeline completed for all chromosomes (including concat/cleanup).\n")
    return per_chr_qc







cfg=read_config("/Users/JJOHN41/Documents/developing_software/postgwas_underdevelopment/postgwas/data/harmonisation/example_input_file.csv")
sample_column_dict=cfg[0]
default_cfg = load_config("harmonisation.yaml")


default_eaf = default_cfg["default_eaf"]
default_external_eaf = default_cfg["default_external_eaf"]
default_info = default_cfg["default_info"]
default_eaf_column = default_cfg["default_eaf_column"]
default_info_column = default_cfg["default_info_column"]
comparison_cols = default_cfg["default_comparison_af_column"]
default_extrnal_eaf = default_cfg["default_comparison_af_file"]
dbsnp = default_cfg["default_dbsnp"]
required_columns_in_sumstat = default_cfg["required_columns_in_sumstat"]
optional_columns_in_sumstat = default_cfg["optional_columns_in_sumstat"]
gwastovcf_main_script_path=default_cfg["gwastovcf_main_script_path"]


# ----------------------------------------------------------------------
# 5. EXAMPLE USAGE
# ----------------------------------------------------------------------
if __name__ == "__main__":
    qc_results = gwas_to_vcf_parallel(
        sumstat_file=sample_column_dict["sumstat_file"],
        sample_column_dict=sample_column_dict,
        output_dir=sample_column_dict['output_folder'],
        resource_folder=sample_column_dict['resourse_folder'],
        default_eaf_file=default_eaf,
        default_info_file=default_info,
        default_comparison_af_file=default_extrnal_eaf,
        dbsnp=dbsnp,
        required_columns_in_sumstat= required_columns_in_sumstat,
        optional_columns_in_sumstat = optional_columns_in_sumstat,
        max_workers = 5,
        user_eaf_file = sample_column_dict['eaffile'],
        user_eaf_column = sample_column_dict['eafcolumn'],
        user_info_file = sample_column_dict['infofile'],
        user_info_column = sample_column_dict['infocolumn'],
        default_eaf_column=default_eaf_column,
        default_comparison_af_column=comparison_cols,
        default_info_column=default_info_column,
        grch37_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh37_check_file.tsv",
        grch38_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh38_check_file.tsv",
        gwastovcf_main_script_path=gwastovcf_main_script_path
    )
    postgwas_qc_summary_df=postgwas_qc(qc_results=qc_results,sample_gwas_dict=sample_column_dict)
    vcf_summary=vcf_qc_summary_polars_full(
        vcf_path=vcf_path,
        tsv_path=tsv_path,
    )

    print("✅ QC summary:")
    for c, qc in qc_results.items():
        print(f"  chr{c}: {list(qc.keys())}")




sumstat_file=sample_column_dict["sumstat_file"]
sample_column_dict=sample_column_dict
output_dir=sample_column_dict['output_folder']
resource_folder=sample_column_dict['resourse_folder']
default_eaf_file=default_eaf
default_info_file=default_info
default_comparison_af_file=default_extrnal_eaf
dbsnp=dbsnp
required_columns_in_sumstat= required_columns_in_sumstat
optional_columns_in_sumstat = optional_columns_in_sumstat
max_workers = 5
user_eaf_file = sample_column_dict['eaffile']
user_eaf_column = sample_column_dict['eafcolumn']
user_info_file = sample_column_dict['infofile']
user_info_column = sample_column_dict['infocolumn']
default_eaf_column=default_eaf_column
default_comparison_af_column=comparison_cols
default_info_column=default_info_column
grch37_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh37_check_file.tsv"
grch38_file = f"{sample_column_dict['resourse_folder']}/GRCh37_38_check_files/GRCh38_check_file.tsv"
gwastovcf_main_script_path=gwastovcf_main_script_path





chromosome=chrom
chr_file=chr_file
sample_column_dict=sample_column_dict.copy()
resource_folder=resource_folder
grch_version=grch_version
user_eaf_file=user_eaf_file
default_eaf_file=default_eaf_file
default_comparison_af_file=default_comparison_af_file
user_info_file=user_info_file
default_info_file=default_info_file
user_eaf_column=user_eaf_column
default_eaf_column=default_eaf_column
default_comparison_af_column=default_comparison_af_column
user_info_column=user_info_column
default_info_column=default_info_column
dbsnp=dbsnp
output_dir=str(output_dir)
gwastovcf_main_script_path=gwastovcf_main_script_path















user_eaf_file='NA'
user_eaf_column = "NA"
user_info_file = "NA"
user_info_column = "NA"
default_eaf_column="EUR"
default_comparison_af_column="INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS"
default_info_column="INFO"

default_eaf="panukb" # ['fingen','panukb','1000G']
default_extrnal_eaf='ALFA'
default_info="panukb" # ['fingen','panukb']
grch_version="GRCh37" # ["GRCh37","GRCh38"]
dbsnp="dbSNP157"

resource_folder="/Users/JJOHN41/Documents/software_resources/resourses/gwas2vcf/"
grch37_file=f"{resource_folder}/GRCh37_38_check_files/GRCh37_check_file.tsv"
grch38_file=f"{resource_folder}/GRCh37_38_check_files/GRCh38_check_file.tsv"

output_dir="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/"



inut_file="filename.csv" # contain following columns sumstat_file	gwas_outputname	chr_col	pos_col	snp_id_col	ea_col	oa_col	eaf_col	beta_or_col	se_col	imp_z_col	pval_col	ncontrol_col	ncase_col	ncontrol	ncase	imp_info_col	delimiter	infofile	infocolumn	eaffile	eafcolumn	liftover	chr_pos_col	resourse_folder	output_folder
default_eaf="panukb" # ['fingen','panukb','1000G']
default_extrnal_eaf='ALFA'
default_info="panukb" # ['fingen','panukb']
grch_version="GRCh37" # ["GRCh37","GRCh38"]
dbsnp="dbSNP157"
required_columns_in_sumstat = [ "CHR", "BP", "A1", "A2", "SNP", "BETA", "SE", "Z", "LP",'FRQ', "N_CON", "NEFF"]
optional_columns_in_sumstat =  ["N_CAS", "INFO"]
output_dir="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/"
user_eaf_file='NA'
user_eaf_column = "NA"
user_info_file = "NA"
user_info_column = "NA"
default_eaf_column="EUR"
default_info_column="INFO"
default_comparison_af_column="INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS"



























def run_pipeline(config_path: str):
    cfg = read_config(config_path)
    df = read_sumstats(cfg["sumstats_path"])

    df = strip_whitespace(df)
    df = split_chr_pos(df)
    validate_gwas_config(df)
    df = fix_chr_column(df)
    df = coerce_numeric(df, ["pos", "p", "eaf", "BETA", "OR"])
    df = check_eaf_or_maf(df)
    df = compute_effective_n(df)

    df = parallel_split_and_run(df, lambda d: d, threads=int(cfg.get("threads", 4)))

    ref_path = load_reference_build(cfg.get("build", "GRCh37"))
    out_tsv = cfg.get("out_tsv", "harmonized.tsv")
    write_tsv(df, out_tsv)

    run_gwas2vcf(out_tsv, out_tsv.replace(".tsv", ""), ref_path, cfg.get("build", "GRCh37"))
    print("✅ Harmonization pipeline completed successfully.")



requiredcolumn_dict = { "chrom": None, "pos": None, "ref": None,"alt": None,
    "b": None,"se": None, "imp_z": None,"nlog_pval": None,
    "alt_freq": None, "rsid": None,
    "n": None,"ncase": None,"ncontrol": None,"neff": None,
    "imp_info": None, 'delimiter':None}


required_columns_in_sumstat = ['chr_col', 'pos_col', 'snp_id_col', 'ea_col', 'oa_col',
    'eaf_col', 'beta_or_col', 'se_col', 'imp_z_col', 'pval_col', 'ncontrol_col','neff_col' ]

optional_columns_in_sumstat = ['ncase_col', 'imp_info_col']


default_eaf="panukb" # ['fingen','panukb','1000G']
default_extrnal_eaf='ALFA'
default_info="panukb" # ['fingen','panukb']
grch_version="GRCh37" # ["GRCh37","GRCh38"]
dbsnp="dbSNP157"

resource_folder="/Users/JJOHN41/Documents/software_resources/resourses/gwas2vcf/"
grch37_file=f"{resource_folder}/GRCh37_38_check_files/GRCh37_check_file.tsv"
grch38_file=f"{resource_folder}/GRCh37_38_check_files/GRCh38_check_file.tsv"

output_dir="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/"



## read_sumstats
cfg=read_config("/Users/JJOHN41/Documents/developing_software/postgwas/data/harmonisation/example_input_file.csv",
                output_dir="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/")

sample_column_dict=cfg[0]

df=read_sumstats(sumstat_file=sample_column_dict['sumstat_file'],
              column_dict=sample_column_dict, 
              output_dir="/Users/JJOHN41/Documents/developing_software/postgwas/data/oudir/")

validate_gwas_config(sample_column_dict, df)

genome_build_info=genome_build(df, grch37_file, grch38_file, sample_column_dict)
grch_version=genome_build_info['inferred_build']

chr_files=split_chr_pos(df, sample_gwas_dict=sample_column_dict, output_folder=output_dir) 



chr="1"

default_eaf_file=f"{resource_folder}/{grch_version}/default_af/tab_files/{grch_version}_{default_eaf}_freq_chr{chr}.tsv.gz"
default_eaf_eafcolumn="EUR"  #[ "AFR","EAS","EUR"]
external_eaf_file='NA' # IF USER PROVIDE HIS/HER OWN SUMMARY STATISTICS EAF FILE

default_info_file=f"{resource_folder}/{grch_version}/default_infoscore/tab_files/{grch_version}_{default_info}_infoscore_chr{chr}.tsv.gz"
default_info_column="INFO"
info_file="NA" # IF USER PROVIDE HIS/HER OWN info FILE
info_column="NA" # IF USER PROVIDE HIS/HER OWN info FILE

genome_fasta_file=f"{resource_folder}/{grch_version}/fasta_files/{grch_version}_chr{chr}.fa" 
default_external_freq_file=f"{resource_folder}/{grch_version}/external_af/vcf_files/{grch_version}_{default_extrnal_eaf}_freq_chr{chr}.vcf.gz"
default_dbsnp_file=f"{resource_folder}/{grch_version}/dbSNP/vcf_files/{grch_version}_{dbsnp}_chr{chr}.vcf.gz" 


chr_file=[x for x in chr_files if f"chr{chr}.tsv" in x][0]

df=pl.read_csv(chr_file, separator="\t")

df=fix_chr_pos_column(df=df, sample_column_dict=sample_column_dict, drop_mt = True)



##identify the default files 
df, eaf_qc_info,sample_column_dict=add_or_calculate_eaf( df=df, sample_column_dict=sample_column_dict,
              eaffile=external_eaf_file,default_eaf_file=default_eaf_file,default_eaf_eafcolumn= default_eaf_eafcolumn)

df,samplesize_qc_info,sample_column_dict=harmonize_sample_sizes( df=df,sample_column_dict=sample_column_dict)

df, qc_info,sample_column_dict=calculate_beta_and_se_from_z(df=df, sample_column_dict=sample_column_dict)

df, qc_info,sample_column_dict=is_beta_or_or(df=df, sample_column_dict=sample_column_dict)

df, qc_info,sample_column_dict=detect_and_convert_pval( df=df, sample_column_dict=sample_column_dict,
    output_col = 'LP',proportion_threshold= 0.10) 

df, qc_info,sample_column_dict=calculate_se_from_beta_pvalue(df=df,sample_column_dict=sample_column_dict)

df, qc_info,sample_column_dict=calculate_z_from_beta_se(df=df, sample_column_dict=sample_column_dict)

df, qc_infoscore,sample_column_dict=add_or_calculate_info( df=df,sample_column_dict=sample_column_dict,
    default_info_file=default_info_file,info_file=info_file,
    default_info_column = default_info_column,info_column= info_column)

df, sample_column_dict=add_or_fix_snp_column(df=df, sample_column_dict=sample_column_dict)


summary_df=export_gwas_sumstat( df=df, sample_column_dict=sample_column_dict,output_dir=output_dir,
    gwas_outputname=sample_column_dict['gwas_outputname'],chr=chr,
    required_columns_in_sumstat=required_columns_in_sumstat,
    optional_columns_in_sumstat=optional_columns_in_sumstat ) 


run_bcftools_munge( output_dir=output_dir, gwas_outputname=sample_column_dict['gwas_outputname'],
    chr=chr,genome_fasta_file=genome_fasta_file )



df_clean, qc_info,sample_column_dict=detect_and_convert_pval( df=df, sample_column_dict=sample_column_dict,
    output_col = sample_column_dict["pval_col"],
    proportion_threshold= 0.10) 



df_clean, qc_info,sample_column_dict=detect_and_convert_pval( df=df, sample_column_dict=sample_column_dict,
    output_col = 'LP',
    proportion_threshold= 0.10) 


