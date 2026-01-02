sample_id="PGC3_SCZ_european"
base_dir='/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"

## convert summstat to vcf and harmonise
docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml

#         --apply-imputation \

docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -it jibinjv/postgwas:1.3 \
  postgwas pipeline --modules flames \
        --apply-filter \
        --heritability \
        --apply-manhattan \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --vcf ${base_dir}/00_harmonised_sumstat/PGC3_SCZ_european_GRCh37_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/ \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
        --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
        --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
        --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
        --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld_clump_population EUR \
        --heritability_tool ldsc \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --heritability_info-min 0.7 \
        --heritability_maf-min 0.01 \
        --maf 0.001 



    --population EUR \
    --corr_method pearson 
    --r2threshold 0.8 \
    --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
    --gwas2vcf_default_config /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml \
    --imputation_tool pred_ld \
    --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \




docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 bash 













docker run --platform=linux/amd64 \
  -e HOME=/tmp \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -it jibinjv/postgwas:1.1 postgwas heritability \
  --sample_id PGC3_SCZ_european \
  --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/heritability/ \
    --heritability_tool ldsc \
    --ldsc_inut /Users/JJOHN41/Documents/developing_software/data/oudir/06_formatter/PGC3_SCZ_european_ldsc_input.tsv \
    --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
    --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
    --samp-prev 0.5 \
    --pop-prev 0.01 \
    --heritability_info-min 0.7 \
    --heritability_maf-min 0.01 



docker run --platform=linux/amd64  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -it jibinjv/postgwas:1.1 postgwas   flames  \
  --finemap_cred_dir /Users/JJOHN41/Documents/developing_software/data/oudir/08_finemap/flames_input/ \
  --magma_genes_out /Users/JJOHN41/Documents/developing_software/data/oudir/09_magma/PGC3_SCZ_european_magma_35up_10down.genes.out \
  --magma_tissue_covar_results /Users/JJOHN41/Documents/developing_software/data/oudir/10_magma_covar/PGC3_SCZ_european.gsa.out \
  --pops_score_file /Users/JJOHN41/Documents/developing_software/data/oudir/11_pops/PGC3_SCZ_european_pops.preds \
  --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
  --sample_id ${sample_id} \
  --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/12_flames_test/




### gcp check 

## convert summstat to vcf and harmonise
docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -e HOME=/tmp \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -it jibinjv/postgwas:1.1 postgwas  harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /mnt/disks/sdd/postgwas_analysis/example_input_file.csv \
        --defaults /mnt/disks/sdd/postgwas_analysis/harmonisation.yaml


sample_id="PGC3_SCZ_european"
base_dir='/mnt/disks/sdd/postgwas_analysis/analysis_results/'
resourse_folder="/mnt/disks/sdd/resourses/postgwas/"
genome_version="GRCh38"


docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -e HOME=/tmp \
  -v /mnt/disks/sdd/:/mnt/disks/sdd/  \
  -it jibinjv/postgwas:1.1 \
  postgwas pipeline \
    --modules flames \
    --heritability \
        --apply-filter \
        --apply-manhattan \
        --apply-imputation \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --vcf ${base_dir}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/ \
        --imputation_tool pred_ld \
        --ref_ld ${resourse_folder}/imputation/pred-ld/ref/ \
        --gwas2vcf_resource ${resourse_folder}/gwas2vcf/ \
        --gwas2vcf_default_config /mnt/disks/sdd/postgwas_analysis/harmonisation.yaml \
        --r2threshold 0.8 \
        --maf 0.001 \
        --population EUR \
        --ref TOP_LD \
        --corr_method pearson \
        --heritability_tool ldsc \
        --merge-alleles ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --heritability_info-min 0.7 \
        --heritability_maf-min 0.01 \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
        --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
        --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
        --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
        --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \





docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.1 postgwas finemap \
  --sample_id test \
  --outdir ${base_dir}/finemap \
  --finemap_method susie \
  --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
  --finemap_input_file ${base_dir}/06_formatter/PGC3_SCZ_european_finemap.tsv \
  --locus_file ${base_dir}/07_ld_clump/PGC3_SCZ_european_LDpruned_EUR_sig.tsv




docker run --platform=linux/amd64   -u $(id -u):$(id -g)   -v /Users/JJOHN41/:/Users/JJOHN41/   -it jibinjv/postgwas:1.3 bash 
eval "$(micromamba shell hook --shell=bash)" 
micromamba activate postgwas  
micromamba activate enricher


python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/finemap/cli.py \
  --sample_id PGC3_SCZ_european \
  --outdir /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/test \
  --locus_file /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/04_ld_clump/PGC3_SCZ_european_LDpruned_EUR_sig.tsv \
  --finemap_method finemap \
  --finemap_ld_ref /Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
  --finemap-in-files /Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/03_formatter/PGC3_SCZ_european_finemap_finemap.tsv --sss


args.outdir="/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/test/"
args.locus_file='/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/04_ld_clump/PGC3_SCZ_european_LDpruned_EUR_sig.tsv'
args.finemap_in_files='/Users/JJOHN41/Documents/developing_software/data/oudir/PGC3_SCZ_european/03_formatter/PGC3_SCZ_european_finemap_finemap.tsv'
args.finemap_ld_ref='/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001'

n_threads=5


import os 

os.system(f" ldstore --in-files {str(ldstore_master)} --read-only-bgen --write-bcor --n-threads {n_threads}")




chr10_18537267_19716878.master


    ], capture_output=True, text=True))