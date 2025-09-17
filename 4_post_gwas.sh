# Step 1: Statistical fine-mapping

PAINTOR -input ${eur_master_file}  -Zhead  Z -LDname ld -in ${eur_file_dir} -out ${eur_result_dir} -mcmc

# Step 2: Statistical fine-mapping based on cross-ancestry GWASs

PAINTOR -input ${cross_ancestry_master_file}  -Zhead  Z_eas,Z_eur -LDname ldeas,ldeur -in ${cross_ancestry_file_dir} -out ${cross_ancestry_result_dir} -mcmc

# Step 3: Credit sets identification based on one causal variant assumption

PAINTOR -input ${eur_master_file} -Zhead  Z -LDname ld -in ${eur_file_dir} -out ${eur_credit_dir} -enumerate 1 -max_causal 1
PAINTOR -input ${cross_ancestry_master_file}  -Zhead  Z_eas,Z_eur -LDname ldeas,ldeur -in ${cross_ancestry_file_dir} -out ${cross_ancestry_credit_dir} -enumerate 1 -max_causal 1

# Step 4: FLAMES

./magma \
 --bfile {PATH_TO_REFERENCE_PANEL_PLINK} \
 --gene-annot {PATH_TO_MAGMA_ANNOT}.genes.annot \
 --pval {PATH_TO_SUMSTATS}.txt ncol=N \
 --gene-model snp-wise=mean \
 --out {DESIRED_ZSCORE_FILENAME}

 ./magma \
--gene-results {DESIRED_ZSCORE_FILENAME}.genes.raw \
--gene-covar {PATH_TO_DOWNLOADED_GTEx_FILE}/gtex_v8_ts_avg_log2TPM.txt \
--out {DESIRED_TISSUE_RELEVANCE_FILENAME}

python pops.py \
--gene_annot_path {PATH_TO_DOWNLOADED_FEATURES}\pops_features_pathway_naive/gene_annot.txt \
--feature_mat_prefix {PATH_TO_DOWNLOADED_FEATURES}\pops_features_pathway_naive/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix {PATH_TO_GENERATED_MAGMA_Z_SCORES}\{DESIRED_ZSCORE_FILENAME} \
--out_prefix {DESIRED_POPS_OUTPUT_PREFIX)

python FLAMES.py annotate \
-o {DESIRED_OUTPUT_DIRECTORY} \
-a {PATH_TO_THE_DOWNLOADED_ANNOTATION_DATA_DIRECTORY} \
-p {DESIRED_POPS_OUTPUT_PREFIX}.preds \
-m {DESIRED_ZSCORE_FILENAME}.genes.out \
-mt {DESIRED_TISSUE_RELEVANCE_FILENAME}.gsa.out \
-id {PATH_TO_INDEXFILE} 

python FLAMES.py FLAMES \
-id {INDEX_FILE_NCLUDING_COLUMN Annotfiles} \
-o {DESIRED_OUTPUT_DIRECTORY}




















