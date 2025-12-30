


# Step 1: Extract Genotypes


# Example using PLINK2 to extract a range and convert to BGEN (preferred by LDstore2)
plink2 --bfile reference_panel \
       --chr 1 --from-bp 1000000 --to-bp 2000000 \
       --export bgen-1.2 bits=8 \
       --out locus_genotypes
       

## Run LDstore2 to compute LD matrix
ldstore --bgen locus_genotypes.bgen \
        --bcor locus_ld.bcor \
        --n-threads 4 \
        --samples samples.txt # Optional
        

# Convert to Matrix Format
ldstore --bcor locus_ld.bcor \
        --merge 4 \
        --matrix locus_ld.matrix
        
        
        
# Prepare FINEMAP Input Files
    # File 1: The Z File (.z)
    # rsid chromosome position allele1 allele2 maf beta se
    # rs123 1 10001 A G 0.2 0.05 0.01
    # rs456 1 10005 T C 0.4 -0.02 0.01

# File 2: The LD File (.ld)
    # z;ld;snp;config;cred;log;n_samples
    # data.z;locus_ld.matrix;data.snp;data.config;data.cred;data.log;50000

# The Master File (.master)# z;ld;snp;config;cred;log



# Run FINEMAP
finemap --sss \
        --in-files finemap_input.files \
        --log finemap.log \
        --n-causal-snps 5 \
        --out finemap_results


# Basic run command
finemap --sss --in-files master_file.master --dataset 1 --log finemap.log --n-causal-snps 5 --out finemap_results