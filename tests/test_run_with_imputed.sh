
## steps
1) harmonisation 
2) sumstat filtering
3) formatter 
4) imputation
5) harmonisation
6) sumstat filtering
7) formatter
8) magma
9) magmacovar
10) pops
11) finemap 
12) flames 


# docker buildx build --no-cache  --platform=linux/amd64 -t jibinjv/postgwas:1.0 .

docker buildx build --platform=linux/amd64 -t jibinjv/postgwas:1.0 .



sample_id="ADHD2022_iPSYCH_deCODE_PGC"
base_dir='/Users/JJOHN41/Documents/developing_software/data/oudir/ADHD2022_iPSYCH_deCODE_PGC/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas harmonisation \
      --nthreads 10 \
      --max-mem 50G \
      --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
      --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas sumstat_filter \
        --nthreads 10 \
        --max-mem 30G \
        --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh38_merged.vcf.gz \
        --sample_id ${sample_id}_GRCh38 \
        --outdir ${base_dir}/2_filtered/ \
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


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas formatter \
    --nthreads 10 \
    --max-mem 30G \
    --seed 10 \
    --vcf ${base_dir}/2_filtered/${sample_id}_GRCh38_filtered.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/3_format_inputs_GRCh38 \
    --format ldpred

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas imputation \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --sample_id ${sample_id} \
        --predld_input_dir ${base_dir}/3_format_inputs_GRCh38/ldpred/ \
        --outdir ${base_dir}/4_imputation_analysis_GRCh38/ \
        --imputation_tool pred_ld \
        --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
        --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
        --r2threshold 0.8 \
        --maf 0.001 \
        --population EUR \
        --ref TOP_LD \
        --corr_method pearson

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas harmonisation \
      --nthreads 10 \
      --max-mem 50G \
      --config ${base_dir}/4_imputation_analysis_GRCh38/${sample_id}_gwas2vcf_config.csv \
      --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas annot_ldblock \
        --nthreads 10 \
        --max-mem 30G \
        --genome-version ${genome_version} \
        --ld_block_population EUR AFR EAS \
        --vcf ${base_dir}/imputed_harmonised/${sample_id}_imputed/1_harmonisation/${sample_id}_imputed_${genome_version}_merged.vcf.gz \
        --ld-region-dir ${resourse_folder}ld_blocks/ 


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas sumstat_filter \
        --nthreads 10 \
        --max-mem 30G \
        --vcf ${base_dir}/imputed_harmonised/${sample_id}_imputed/1_harmonisation/${sample_id}_imputed_${genome_version}_merged.vcf.gz \
        --sample_id ${sample_id}_GRCh37 \
        --outdir ${base_dir}/2_filtered_imputed/ \
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


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas formatter \
    --nthreads 10 \
    --max-mem 30G \
    --seed 10 \
    --vcf ${base_dir}/2_filtered_imputed/${sample_id}_${genome_version}_filtered.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/3_format_inputs \
    --format magma finemap ldpred ldsc


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas heritability \
        --nthreads 4 \
        --max-mem 30G \
        --seed 100 \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/4_ldsc_analysis/ \
        --sumstats ${base_dir}/3_format_inputs/${sample_id}_ldsc_input.tsv \
        --heritability_tool ldsc \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --info-min 0.7 \
        --maf-min 0.01 \
        --docker-image jibinjv/ldsc:1.0.1 \
        --platform linux/amd64

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas ld_clump \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --sample_id ${sample_id} \
    --vcf ${base_dir}/2_filtered_imputed/${sample_id}_${genome_version}_filtered.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/5_ldclump_analysis/ \
    --ld-mod by_regions \
    --population EUR 


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas magma \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --genome-version ${genome_version} \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/7_magma_gene_pathway_assoc/ \
    --snp_loc_file ${base_dir}/3_format_inputs/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/3_format_inputs/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
    --num_batches 6 \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas magmacovar \
      --nthreads 6 \
      --max-mem 30G \
      --seed 10 \
      --sample_id ${sample_id} \
      --outdir ${base_dir}/7_magma_gene_pathway_assoc/ \
      --covariates ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
      --covar_model condition-hide=Average \
      --covar_direction greater \
      --magama_gene_assoc_raw ${base_dir}/7_magma_gene_pathway_assoc/${sample_id}_magma_35up_10down.genes.raw 


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas pops \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/8_pops_analysis/ \
    --magma_assoc_prefix ${base_dir}/7_magma_gene_pathway_assoc/${sample_id}_magma_35up_10down \
    --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
    --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas finemap \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/9_finemap/ \
    --locus_file ${base_dir}/5_ldclump_analysis/${sample_id}_LDpruned_EUR_sig.tsv \
    --finemap_input_file ${base_dir}/3_format_inputs/${sample_id}_finemap.tsv \
    --finemap_method susie \
    --lp_threshold 7.3 \
    --L 10 \
    --min_ram_per_worker_gb 10 \
    --timeout_ld_seconds 300 \
    --timeout_susie_seconds 300 \
    --finemap_skip_mhc \
    --finemap_mhc_start 25000000 \
    --finemap_mhc_end 35000000 \
    --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas flames \
  --pops ${base_dir}/8_pops_analysis/${sample_id}_pops.preds \
  --magma_tissue_covar_results ${base_dir}/7_magma_gene_pathway_assoc/${sample_id}.gsa.out \
  --finemap_cred_dir ${base_dir}/9_finemap/falmes_inputput/ \
  --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
  --magma_genes_out ${base_dir}/7_magma_gene_pathway_assoc/${sample_id}_magma_35up_10down.genes.out \
  --out ${base_dir}/10_flames_analysis



docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas manhattan \
      --nthreads 10 \
      --max-mem 30G \
      --seed 10 \
      --vcf ${base_dir}/2_filtered_imputed/${sample_id}_${genome_version}_filtered.vcf.gz \
      --sample_id ${sample_id} \
      --outdir ${base_dir}/11_manhatton \
      --pdf ${base_dir}/11_manhatton/${sample_id}_manhattonplot.pdf \
      --genome GRCh37 \
      --pheno ${sample_id}_imputed \
      --csq \
      --nauto 22 \
      --min-af 0.001 \
      --min-lp 1 \
      --loglog-pval 10 \
      --cyto-ratio 25 \
      --spacing 8 \
      --width 10.0 \
      --fontsize 10



























































































































































docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas heritability \
        --nthreads 4 \
        --max-mem 30G \
        --seed 100 \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/4_ldsc_analysis/ \
        --sumstats ${base_dir}/3_format_inputs_GRCh38/${sample_id}_ldsc_input.tsv \
        --heritability_tool ldsc \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --info-min 0.7 \
        --maf-min 0.01 \
        --docker-image jibinjv/ldsc:1.0.1 \
        --platform linux/amd64


docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
    postgwas formatter \
    --nthreads 10 \
    --max-mem 30G \
    --seed 10 \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_GRCh38_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/3_format_inputs_GRCh38 \
    --format magma finemap ldpred ldsc

docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents/ \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -it jibinjv/postgwas:1.0 \
  postgwas harmonisation \
      --nthreads 10 \
      --max-mem 50G \
      --config /Users/JJOHN41/Documents/developing_software/data/oudir/ADHD2022_iPSYCH_deCODE_PGC/5_imputation_analysis/ADHD2022_iPSYCH_deCODE_PGC_gwas2vcf_config.csv \
      --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml


postgwas ld_clump \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --sample_id ${sample_id} \
    --bcftools /Users/JJOHN41/bin/bcftools \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_${genome_version}_merged.vcf.gz  \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/6_ldclump_analysis/ \
    --ld-mod by_regions \
    --population EUR 

postgwas magma \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \
    --genome-version ${genome_version} \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/7_magma_gene_pathway_assoc/ \
    --snp_loc_file ${base_dir}/3_format_inputs/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/3_format_inputs/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file  ${resourse_folder}/magma/gene_loc/NCBI37.3/NCBI37.3.gene_withGeneSymbol.loc \
    --geneset_file ${resourse_folder}/msigdb_v2025.1.Hs_GMTs/msigdb.v2025.1.Hs.symbols.gmt \
    --num_batches 6 \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL


postgwas magma \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \
    --genome-version ${genome_version} \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/7_magma_gene_pathway_assoc/ \
    --snp_loc_file ${base_dir}/3_format_inputs/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/3_format_inputs/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file  ${resourse_folder}/magma/gene_loc/NCBI37.3/NCBI37.3.gene_withGeneSymbol.loc \
    --geneset_file ${resourse_folder}/msigdb_v2025.1.Hs_GMTs/msigdb.v2025.1.Hs.symbols.gmt \
    --num_batches 6 \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL


postgwas magma \
    --nthreads 10 \
    --max-mem 60G \
    --seed 10 \
    --magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \
    --genome-version ${genome_version} \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/8_tissue_specificity/ \
    --snp_loc_file ${base_dir}/3_format_inputs/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/3_format_inputs/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
    --num_batches 6 \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL

postgwas magmacovar \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/8_tissue_specificity/ \
    --magma /Users/JJOHN41/Documents/software_resources/softwares/magma_v1.10_mac/magma \
    --covariates ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
    --covar_model condition-hide=Average \
    --covar_direction greater \
    --magama_gene_assoc_raw ${base_dir}/8_tissue_specificity/${sample_id}_magma_35up_10down.genes.raw 

postgwas pops \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/8_pops_analysis/ \
    --magma_assoc_prefix ${base_dir}/8_tissue_specificity/${sample_id}_magma_35up_10down \
    --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
    --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \

postgwas finemap \
    --nthreads 6 \
    --max-mem 30G \
    --seed 10 \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/9_finemap/ \
    --locus_file ${base_dir}/6_ldclump_analysis/${sample_id}_LDpruned_EUR_sig.tsv \
    --finemap_input_file ${base_dir}/3_format_inputs/${sample_id}_finemap.tsv \
    --finemap_method susie \
    --lp_threshold 7.3 \
    --L 10 \
    --min_ram_per_worker_gb 10 \
    --timeout_ld_seconds 300 \
    --timeout_susie_seconds 300 \
    --finemap_skip_mhc \
    --finemap_mhc_start 25000000 \
    --finemap_mhc_end 35000000 \
    --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --plink /Users/JJOHN41/Documents/software_resources/softwares/plink 


postgwas flames \
 --pops ${base_dir}/8_pops_analysis/${sample_id}_pops.preds \
 --magma_tissue_covar_results ${base_dir}/8_tissue_specificity/${sample_id}.gsa.out \
 --finemap_cred_dir ${base_dir}/9_finemap/falmes_inputput/ \
 --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
 --magma_genes_out ${base_dir}/8_tissue_specificity/${sample_id}_magma_35up_10down.genes.out \
 --out ${base_dir}/10_flames_analysis



postgwas manhattan \
    --nthreads 10 \
    --max-mem 30G \
    --seed 10 \
    --vcf ${base_dir}/1_harmonisation/${sample_id}_${genome_version}_merged.vcf.gz \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/11_manhatton \
    --pdf ${base_dir}/11_manhatton/${sample_id}_manhattonplot.pdf \
    --genome GRCh37 \
    --pheno ${sample_id} \
    --csq \
    --nauto 22 \
    --min-af 0.001 \
    --min-lp 1 \
    --loglog-pval 10 \
    --cyto-ratio 25 \
    --spacing 8 \
    --width 10.0 \
    --fontsize 10
