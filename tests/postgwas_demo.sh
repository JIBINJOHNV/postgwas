

sample_id="PGC3_SCZ_european"
base_dir='/Users/JJOHN41/Documents/developing_software/data/oudir/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh38"

## convert summstat to vcf and harmonise
docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -it jibinjv/postgwas:1.0 postgwas  harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml


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
        --vcf ${base_dir}/00_harmonised_sumstat/${sample_id}_GRCh38_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/ \
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
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001






folder_path='/Users/JJOHN41/Documents/developing_software/data/oudir/03_imputation/chrwise/'
folder_path = Path(folder_path)


from postgwas.imputation.pred_ld.pred_ld_runner import (
             run_pred_ld_parallel,process_pred_ld_results_all_parallel         )


process_pred_ld_results_all_parallel(
    folder_path=folder_path,
    output_path=output_path,
    gwas2vcf_resource_folder=gwas2vcf_resource_folder,
    output_prefix=output_prefix,
    corr_method= "pearson",
    threads= 1,
)


    --heritability --apply-filter --apply-imputation --apply-manhattan 

  /Users/JJOHN41/Documents/developing_software/data/oudir/03_imputation/chrwise
