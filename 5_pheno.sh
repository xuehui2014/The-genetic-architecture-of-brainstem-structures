# Step 1: Genetic correlation analyses

ldsc.py --rg ${gwas1}.sumstats.gz,${gwas2}.sumstats.gz  --ref-ld-chr ${ld_file} --w-ld-chr ${ld_file} --out ${gwas1}_${gwas2}

# Step 2: Genetic colocalization analyses

Rscript coloc_code.R

# Step 3: conjFDR analyses

matlab -nodisplay -nosplash < runme_gwas1_gwas2.m


