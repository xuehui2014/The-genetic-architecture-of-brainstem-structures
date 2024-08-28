########## Step1: Brainstem segmentation ##########

subj=010001
segmentBS.sh  ${subj}

########## Step2: Genetic data QC ##########

##########a. Sample-level QC in ABCD ##########

python inferancestry.py --par ./abcd.par  ###keep participants with propotion of European genetic ancestry > 80%

########## b. ABCD_PCA ##########

Rscript ABCD_GRM.r --bfile ABCD --force

######### c. marker-level QC ##########

genetic_file={CHIMGEN,ABCD,UKBB}
plink --bfile ${genetic_file} --maf 0.005  --write-snplist  -snps-only just-acgt --out maf_snp ### maf > 0.005 and exclude non-SNP
plink2  --bfile  ${genetic_file} --export bgen-1.3   --out ${genetic_file}
qctool  -g ${genetic_file}.bgen  -snp-stats -osnp   ${genetic_file}_hwe_qctool ### PHWE ≥ 1 × 10-7

######### d. make genetic file  ##########

plink --bfile ${genetic_file} --keep qualified_sample --extract qualified_snplist --make-bed --out "${genetic_file}_qualified" ##extract MAF ≥ 0.5%, imputation quality r2 > 0.3 and PHWE ≥ 1 × 10-7 and keep qualified participents afted sample-level QC

######### Step3: Preperation for GWAS ##########

######### a. Combat harmonization  ##########

a1=load('input_trait')
b1=load('input_center_number')
c1=load('cov_sex_age')
harmonization_run(a1,b1,c1);

######### a. make grm file  ##########

for i in {1..50}
do
gcta  --mbfile genetic_file_list --keep qualified_sample --make-grm-part 50 $i --thread-num 50 --out genetic_grm 
done

cat genetic_grm.*.grm.id > genetic.grm.id
cat genetic_grm.*.grm.bin > genetic.grm.bin
cat genetic_grm.*.grm.N.bin > genetic.grm.N.bin

gcta --grm genetic --make-bK-sparse 0.05 --out genetic_sp_grm


######### b. quantile normalization  ##########

plink --bfile ${genetic_file} --keep qualified_sample --extract qualified_snplist --pheno input_${phenotype_file} --covar ${cov_file} --quantile-normalize --write-covar --out ${phenotype_file}

mv ${phenotype_file}.cov ${phenotype_file}

######### Step4: GWAS ##########

######### a. Whole brainstem volume and substructure abbsolute volume  ##########
for i in {1..5}
do

gcta --bfile ${genetic_file}    \
 --fastGWA-mlm  --pheno ${phenotype_file}  \
--qcovar ${qcov_file_no_brainstem}   --covar ${cov_file}    --maf 0.005 --mpheno $i \
 --grm-sparse  genetic_sp_grm \
--threads 50           \
--out result_included_$i


gcta --bfile ${genetic_file}    \
 --fastGWA-mlm  --pheno ${phenotype_file}  \
--model-only \
--qcovar ${qcov_file_no_brainstem}        --covar ${cov_file}    --maf 0.005  --mpheno $i   \
 --grm-sparse  genetic_sp_grm     --threads 50      \
--out result_included_$i


gcta --bfile ${X_genetic_file}   --chr 23 \
--pheno ${phenotype_file} --covar ${cov_file} \
--qcovar ${qcov_file_no_brainstem}       --maf 0.005       --mpheno $i     --threads 50             \
--load-model "result_included_$i.fastGWA"   \
--out X_result_included_$i

done

######### b. Substructure relative volume ##########

for i in {2..5}
do

gcta --bfile ${genetic_file}    \
 --fastGWA-mlm  --pheno ${phenotype_file}  \
--qcovar ${qcov_file_brainstem}    --covar ${cov_file}    --maf 0.005 --mpheno $i \
 --grm-sparse  genetic_sp_grm \
--threads 50           \
--out result_excluded_$i


gcta --bfile ${genetic_file}    \
 --fastGWA-mlm  --pheno ${phenotype_file}  \
--model-only \
--qcovar ${qcov_file_brainstem}         --covar ${cov_file}    --maf 0.005  --mpheno $i   \
 --grm-sparse  genetic_sp_grm     --threads 50      \
--out result_excluded_$i


gcta --bfile ${X_genetic_file}   --chr 23 \
--pheno ${phenotype_file} --covar ${cov_file} \
--qcovar ${qcov_file_brainstem}        --maf 0.005       --mpheno $i     --threads 50             \
--load-model "result_excluded_$i.fastGWA"   \
--out X_result_excluded_$i

######### c. Within-EUR meta-analyses ##########

metal meta_eur.txt

######### d. Cross-ancestry meta-analyses ##########

metal meta_cross_ancestry.txt

######### Step5: Post-GWAS analyses ##########

######### a. pooling #########

plink --bfile ${genetic_file}   --clump ${gwas_file}             --clump-p1 5.56e-9            --clump-p2 5e-8       --clump-r2 0.1         --clump-kb 500         --clump-field "P"   --out p_5e-8_9_r2_0.1_window_500_${gwas_file}

######### b. h2 #########

######### I. cross-ancestry LD matrix  #########

plink --bfile ${cross_ancestry_genetic_file} --exclude high-LD-regions.txt --range --out ${cross_ancestry_genetic_file_rm_high-LD}
plink --bfile ${cross_ancestry_genetic_file_rm_high-LD} --indep-pairwise 50 5 0.2 --out file
plink --bfile ${cross_ancestry_genetic_file_rm_high-LD} --extract file.prune.in --pca 40 --out cross_ancestry_pca
python cov-ldsc/ldsc.py  --bfile   ${cross_ancestry_genetic_file} --cov cross_ancestry_pca.eigenvec --l2   --ld-wind-cm 20  --out ${cross_ancestry_genetic_file_ld}

######### II. h2 calculation  #########

munge_sumstats.py --sumstats ${gwas_file} --out ${gwas_file} --merge-alleles w_hm3.snplist --chunksize 500000
ldsc.py --h2 ${gwas_file}.sumstats.gz --ref-ld-chr ${ld_file} --w-ld-chr ${ld_file} --out  ${gwas_file}

######### c. Cross-ancestry genetic correlation #########

popcorn fit -v 1 --cfile ${scores_file} --gen_effect --sfile1 ${eas_gwas_file} --sfile2 ${eur_gwas_file} "output/${eas_gwas_file}_${eur_gwas_file}_eff"

popcorn fit -v 1 --cfile ${scores_file}  --gen_effect --sfile1 ${eas_gwas_file} --sfile2 ${eur_gwas_file} "output/${eas_gwas_file}_${eur_gwas_file}_imp"

######### d. Statistical fine-mapping #########

######### I. Statistical fine-mapping based on EUR-GWASs #########

PAINTOR -input ${eur_master_file}  -Zhead  Z -LDname ld -in ${eur_file_dir} -out ${eur_result_dir} -mcmc

######### II. Statistical fine-mapping based on cross-ancestry GWASs #########

PAINTOR -input ${cross_ancestry_master_file}  -Zhead  Z_eas,Z_eur -LDname ldeas,ldeur -in ${cross_ancestry_file_dir} -out ${cross_ancestry_result_dir} -mcmc

######### III. Credit sets  #########

PAINTOR -input ${eur_master_file} -Zhead  Z -LDname ld -in ${eur_file_dir} -out ${eur_credit_dir} -enumerate 1 -max_causal 1
PAINTOR -input ${cross_ancestry_master_file}  -Zhead  Z_eas,Z_eur -LDname ldeas,ldeur -in ${cross_ancestry_file_dir} -out ${cross_ancestry_credit_dir} -enumerate 1 -max_causal 1

######### e. Genes #########

######### I. MAGMA #########

magma   --annotate window=0,0  	--snp-loc 1000G.EUR.bim  	--gene-loc NCBI37.3.gene.loc   	--out step1_EUR
magma     --bfile 1000G.EUR     --pval ${gwas_file} use='SNP,P' N=65621     --gene-annot step1_EUR.genes.annot     --out step2_${gwas_file}

######### II. S-MultiXscan #########

python gwas_parsing.py \
-gwas_file ${gwas_file}.txt.gz \
-liftover hg19ToHg38.over.chain.gz \
-snp_reference_metadata variant_metadata.txt.gz METADATA \
-output_column_map SNP variant_id \
-output_column_map A2 non_effect_allele \
-output_column_map A1 effect_allele \
-output_column_map BETA effect_size \
-output_column_map P pvalue \
-output_column_map CHR chromosome \
--chromosome_format \
-output_column_map POS position \
-output_column_map Freq1 frequency \
--insert_value sample_size 65621 --insert_value n_cases 0 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size \
-output harmonized_${gwas_file}.txt.gz

for CHR in $(seq 1 22);do
for batch in $(seq 0 9);do
python gwas_summary_imputation.py \
-by_region_file eur_ld.bed.gz \
-gwas_file harmonized_${gwas_file}.txt.gz \
-parquet_genotype chr"${CHR}".variants.parquet \
-parquet_genotype_metadata variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome ${CHR} \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch ${batch} \
--standardise_dosages \
-output "${imputation_dir}/${gwas_file}_chr${CHR}_sb${batch}_reg0.1_ff0.01_by_region.txt.gz"
done
done

python gwas_summary_imputation_postprocess1.py \
-gwas_file harmonized_${gwas_file}.txt.gz \
-folder ${imputation_dir}/ \
-pattern ${gwas_file}_chr* \
-parsimony 7 \
-output imputed_${gwas_file}.txt.gz

python SPrediXcan.py \
--gwas_file  imputed_${gwas_file}.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path mashr_${region}.db \
--covariance mashr_${region}.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file ${metaxcan_dir}/${gwas_file}_${region}.csv

python SMulTiXcan.py \
--models_folder ${models_dir} \
--models_name_pattern "mashr_Brain_(.*).db" \
--snp_covariance gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder ${metaxcan_dir} \
--metaxcan_filter "${gwas_file}_Brain_(.*).csv" \
--metaxcan_file_name_parse_pattern "${gwas_file}_Brain_(.*).csv" \
--gwas_file ${gwas_file}.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output ${gwas_file}_smultixcan.csv

######### III. H-MAGMA #########

annotation_file=iPSC_derived_neuro
magma     --bfile g1000_eur     --pval ../h2/input_eur_${gwas} use='SNP,P' N=65621    --gene-annot ${annotation_file}.genes.annot  --out ${gwas}_${annotation_file}

######### f. Shared genetic architectures with other non-imaging phenotypes #########

######### I. Genetic correlation analyses #########

ldsc.py --rg ${gwas1}.sumstats.gz,${gwas2}.sumstats.gz  --ref-ld-chr ${ld_file} --w-ld-chr ${ld_file} --out ${gwas1}_${gwas2}

######### II. Genetic colocalization analyses #########

Rscript coloc_code.R

######### III. conjFDR analyses #########

matlab -nodisplay -nosplash < runme_gwas1_gwas2.m

######### Software and algorithms websites #########

######### 1. Freesurfer brainstem substructures segmentation: https://surfer.nmr.mgh.harvard.edu/fswiki/BrainstemSubstructures
######### 2. SNPWEIGHTS: https://www.hsph.harvard.edu/alkes-price/software/
######### 3. ABCD_GeneticPCs: https://github.com/robloughnan/ABCD_GeneticPCs_and_Relatedness
######### 4. PLINK: https://www.cog-genomics.org/plink/1.9/
######### 5. PLINK2: https://www.cog-genomics.org/plink/2.0/
######### 6. QCTOOL: https://www.chg.ox.ac.uk/~gav/qctool_v2/
######### 7. fastGWA: https://yanglab.westlake.edu.cn/software/gcta/#Overview
######### 8. LDSC: https://github.com/bulik/ldsc
######### 9. cov-LDSC: https://github.com/yang-luo-lab/cov-ldsc
######### 10. Popcorn: https://github.com/brielin/Popcorn
######### 11. PAINTOR: https://github.com/gkichaev/PAINTOR_V3.0
######### 12. MAGMA: https://cncr.nl/research/magma/
######### 13. S-MultiXcan: https://github.com/hakyimlab/MetaXcan
######### 14. H-MAGMA: https://github.com/thewonlab/H-MAGMA
######### 15. coloc: https://github.com/chr1swallace/coloc
######### 16. conjFDR: https://github.com/precimed/pleiofdr
