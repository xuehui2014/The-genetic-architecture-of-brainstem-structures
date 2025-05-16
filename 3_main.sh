# Step 1: GWAS for whole brainstem volume and brainstem substructure absolute volumes
for i in {1..5}
    do
    gcta --bfile ${genetic_file}    \
        --fastGWA-mlm  --pheno ${phenotype_file}  \
        --qcovar ${qcov_file_no_brainstem}   --covar ${cov_file}    --maf 0.005 --mpheno $i \
        --grm-sparse  genetic_sp_grm \
        --out result_included_$i
    gcta --bfile ${genetic_file}    \
        --fastGWA-mlm  --pheno ${phenotype_file}  \
        --model-only \
        --qcovar ${qcov_file_no_brainstem}        --covar ${cov_file}    --maf 0.005  --mpheno $i   \
        --grm-sparse  genetic_sp_grm       \
        --out result_included_$i
    gcta --bfile ${X_genetic_file}   --chr 23 \
        --pheno ${phenotype_file} --covar ${cov_file} \
        --qcovar ${qcov_file_no_brainstem}       --maf 0.005       --mpheno $i                \
        --load-model "result_included_$i.fastGWA"   \
        --out X_result_included_$i
done

# Step 2: GWAS for whole brainstem substructure relative volumes
for i in {2..5}
    do
    gcta --bfile ${genetic_file}    \
        --fastGWA-mlm  --pheno ${phenotype_file}  \
        --qcovar ${qcov_file_brainstem}    --covar ${cov_file}    --maf 0.005 --mpheno $i \
        --grm-sparse  genetic_sp_grm \
        --out result_excluded_$i
    gcta --bfile ${genetic_file}    \
        --fastGWA-mlm  --pheno ${phenotype_file}  \
        --model-only \
        --qcovar ${qcov_file_brainstem}         --covar ${cov_file}    --maf 0.005  --mpheno $i   \
        --grm-sparse  genetic_sp_grm        \
        --out result_excluded_$i
    gcta --bfile ${X_genetic_file}   --chr 23 \
        --pheno ${phenotype_file} --covar ${cov_file} \
        --qcovar ${qcov_file_brainstem}        --maf 0.005       --mpheno $i           \
        --load-model "result_excluded_$i.fastGWA"   \
        --out X_result_excluded_$i
done

# Step 3: GWAS meta analyses to generate EUR-GWAS

metal meta_eur.txt

# Step 4: GWAS meta analyses to generate EUR-EAS-GWAS

metal meta_cross_ancestry.txt

# Step 5: GWAS meta analyses to generate admixed-GWAS

metal meta_admixed_ancestry.txt






        
