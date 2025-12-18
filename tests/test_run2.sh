python src/postgwas/imputation/cli.py \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --sample_id ${sample_id} \
    --predld_input_dir ${base_dir}/3_format_inputs_GRCh38/ldpred/ \
    --outdir ${base_dir}/5_imputation_analysis/ \
    --imputation_tool pred_ld \
    --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
    --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --r2threshold 0.8 \
    --maf 0.001 \
    --population EUR \
    --ref TOP_LD \
    --corr_method pearson



"[bold green]Default:[/bold green] [cyan]0.8[/cyan]"
"[bold bright_red]Required[/bold bright_red]: 


annot_ldblock,sumstat_filter,ld_clump,finemap,magma,magmacovar,pops,flames 



get_magma_binary_parser,
get_bcftools_binary_parser,
get_plink_binary_parser

sumstat_filter , imputation not showing help 


genome_version=GRCh37

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/harmonisation/cli.py \
      --nthreads 10 \
      --max-mem 50G \
      --config /Users/JJOHN41/Documents/developing_software/data/oudir/ADHD2022_iPSYCH_deCODE_PGC/5_imputation_analysis/ADHD2022_iPSYCH_deCODE_PGC_gwas2vcf_config.csv \
      --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/annot_ldblock/cli.py \
    --nthreads 10 \
        --max-mem 30G \
        --genome-version ${genome_version} \
        --ld_block_population EUR AFR EAS \
        --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --sample_id ${sample_id}_${genome_version} \
        --outdir ${base_dir}/2_ldblock_annotated/

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/sumstat_filter/cli.py \
        --nthreads 10 \
        --max-mem 30G \
        --vcf ${base_dir}/2_ldblock_annotated/${sample_id}_${genome_version}_ldblock.vcf.gz \
        --sample_id ${sample_id}_${genome_version} \
        --outdir ${base_dir}/3_filtered/ \
        --pval-cutoff 0.0001 \
        --maf-cutoff 0.001 \
        --allelefreq-diff-cutoff 0.2 \
        --info-cutoff 0.7 \
        --external-af-name EUR \
        --include-indels \
        --include-palindromic \
        --palindromic-af-lower 0.4 \
        --palindromic-af-upper 0.6 \
        --remove-mhc  \
        --mhc-chrom 6 \
        --mhc-start 25000000 \
        --mhc-end 34000000  

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/formatter/cli.py \
    --nthreads 10 \
    --max-mem 30G \
    --seed 10 \
    --vcf ${base_dir}/3_filtered/${sample_id}_${genome_version}_filtered.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/4_downstream_inputs \
    --format magma finemap ldpred ldsc

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/ld_clump/cli.py \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --sample_id ${sample_id} \
    --vcf  ${base_dir}/3_filtered/${sample_id}_${genome_version}_filtered.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/5_ldclump_analysis/ \
    --ld-mod by_regions \
    --population EUR 

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/gene_assoc/cli.py \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --genome-version ${genome_version} \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/6_magma_gene_pathway_assoc/ \
    --snp_loc_file ${base_dir}/4_downstream_inputs/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/4_downstream_inputs/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
    --num_batches 6 \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/magmacovar/cli.py \
      --nthreads 6 \
      --max-mem 30G \
      --seed 10 \
      --sample_id ${sample_id} \
      --outdir ${base_dir}/6_magma_gene_pathway_assoc/ \
      --covariates ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
      --magama_gene_assoc_raw ${base_dir}/6_magma_gene_pathway_assoc/${sample_id}_magma_35up_10down.genes.raw \

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pops/cli.py \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/7_pops_analysis/ \
    --magma_assoc_prefix ${base_dir}/6_magma_gene_pathway_assoc/${sample_id}_magma_35up_10down \
    --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
    --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt 






direction-covar=greater condition-hide=Average



if in pipeline flames module


steps is as followis 

annot_ldblock
sumstat_filter
formatter
ld_clump
magma
magmacovar
pops
flames
