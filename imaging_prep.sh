# Step 1: Brainstem segmentation

subj=010001
segmentBS.sh  ${subj}

# Step 2: Combat harmonization

a1=load('input_trait')
b1=load('input_center_number')
c1=load('cov_sex_age')
harmonization_run(a1,b1,c1);

# Step 3: Quantile normalization

plink --bfile ${genetic_file} --keep qualified_sample --extract qualified_snplist --pheno input_${phenotype_file} --covar ${cov_file} --quantile-normalize --write-covar --out ${phenotype_file}

mv ${phenotype_file}.cov ${phenotype_file}
