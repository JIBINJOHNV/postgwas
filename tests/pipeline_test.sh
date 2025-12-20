
sample_id="ADHD2022_iPSYCH_deCODE_PGC"
base_dir='/Users/JJOHN41/Documents/developing_software/data/oudir/ADHD2022_iPSYCH_deCODE_PGC/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/harmonisation/cli.py \
      --nthreads 10 \
      --max-mem 50G \
      --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
      --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
  --modules annot_ldblock \
  --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
  --genome-version GRCh37 \
  --sample_id ${sample_id} \
  --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing \
  --ld_block_population EUR AFR EAS \
  --ld-region-dir ${resourse_folder}ld_blocks/ 

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
    --modules formatter \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing \
    --format magma finemap ldpred ldsc 


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
    --modules sumstat_filter \
    --genome-version GRCh37 \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
    --modules formatter  --apply-filter  \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing \
    --format magma finemap ldpred ldsc 

python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
    --modules formatter \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh37_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing \
    --format ldpred 


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/pipeline/cli.py \
    --modules finemap \
    --heritability \
        --apply-filter \
        --apply-imputation \
        --apply-manhattan \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --vcf  ${base_dir}/1_harmonisation/${sample_id}_GRCh38_chr1_chr2.vcf.gz \
        --sample_id ${sample_id} \
        --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing/ \
        --imputation_tool pred_ld \
        --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
        --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
        --gwas2vcf_default_config /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml \
        --r2threshold 0.8 \
        --maf 0.001 \
        --population EUR \
        --ref TOP_LD \
        --corr_method pearson \
        --heritability_tool ldsc \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --info-min 0.7 \
        --maf-min 0.01 \
        --docker-image jibinjv/ldsc:1.0.1 \
        --platform linux/amd64 \
        --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
        --magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \
        --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
        --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
        --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --plink /Users/JJOHN41/Documents/software_resources/softwares/plink 


docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -it jibinjv/postgwas:1.0 postgwas  harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml



sample_id="PGC3_SCZ_european"

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 postgwas pipeline \
    --modules flames \
    --heritability \
        --apply-filter \
        --apply-imputation \
        --apply-manhattan \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --vcf /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_EUR/harmonised_susmstat/${sample_id}_GRCh38_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_EUR/ \
        --imputation_tool pred_ld \
        --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
        --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
        --gwas2vcf_default_config /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml \
        --r2threshold 0.8 \
        --maf 0.001 \
        --population EUR \
        --ref TOP_LD \
        --corr_method pearson \
        --heritability_tool ldsc \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --info-min 0.7 \
        --maf-min 0.01 \
        --docker-image jibinjv/ldsc:1.0.1 \
        --platform linux/amd64 \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
        --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
        --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
        --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
        --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \








docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -it jibinjv/postgwas:1.0 postgwas  harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_EUR/03_imputation/PGC3_SCZ_european_gwas2vcf_config.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml
















--plink /Users/JJOHN41/Documents/software_resources/softwares/plink \
--magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \

postgwas pipeline --module sumstat_filter --apply-imputation \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh38_chr1_chr2.vcf.gz \
    --sample_id ${sample_id} \
    --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/pipeline_testing/ \
    --imputation_tool pred_ld \
    --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
    --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --gwas2vcf_default_config /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml \
    --r2threshold 0.8 \
    --maf 0.001 \
    --population EUR 