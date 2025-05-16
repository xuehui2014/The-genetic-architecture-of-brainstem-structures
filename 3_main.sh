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

# Step 6: Pooling-clumping

plink --bfile ${genetic_file}   --clump ${gwas_file}             --clump-p1 5.56e-9            --clump-p2 5e-8       --clump-r2 0.1         --clump-kb 500         --clump-field "P"   --out p_5e-8_9_r2_0.1_window_500_${gwas_file}

# Step 6: Cross-ancestry LD score

plink --bfile ${cross_ancestry_genetic_file} --exclude high-LD-regions.txt --range --out ${cross_ancestry_genetic_file_rm_high-LD}
plink --bfile ${cross_ancestry_genetic_file_rm_high-LD} --indep-pairwise 50 5 0.2 --out file
plink --bfile ${cross_ancestry_genetic_file_rm_high-LD} --extract file.prune.in --pca 40 --out cross_ancestry_pca
python cov-ldsc/ldsc.py  --bfile   ${cross_ancestry_genetic_file} --cov cross_ancestry_pca.eigenvec --l2   --ld-wind-cm 20  --out ${cross_ancestry_genetic_file_ld}


# Step 7: h2 calculation

munge_sumstats.py --sumstats ${gwas_file} --out ${gwas_file} --merge-alleles w_hm3.snplist --chunksize 500000
ldsc.py --h2 ${gwas_file}.sumstats.gz --ref-ld-chr ${ld_file} --w-ld-chr ${ld_file} --out  ${gwas_file}

# Step 8: Cross-ancestry genetic correlation

popcorn fit -v 1 --cfile ${scores_file} --gen_effect --sfile1 ${eas_gwas_file} --sfile2 ${eur_gwas_file} "output/${eas_gwas_file}_${eur_gwas_file}_eff"

popcorn fit -v 1 --cfile ${scores_file}  --gen_effect --sfile1 ${eas_gwas_file} --sfile2 ${eur_gwas_file} "output/${eas_gwas_file}_${eur_gwas_file}_imp"






        
