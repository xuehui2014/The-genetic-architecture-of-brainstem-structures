# Step 1: SNPWeight
python inferancestry.py --par ./abcd.par  
# keep participants with propotion of European genetic ancestry > 80%


# Step 2: PCA calculation
Rscript ABCD_GRM.r --bfile ABCD --force


# Step 3: Marker-level QC
genetic_file={CHIMGEN,ABCD,UKBB}
plink --bfile ${genetic_file} --maf 0.005  --write-snplist  -snps-only just-acgt --out maf_snp ### maf > 0.005 and exclude non-SNP
plink2  --bfile  ${genetic_file} --export bgen-1.3   --out ${genetic_file}
qctool  -g ${genetic_file}.bgen  -snp-stats -osnp   ${genetic_file}_hwe_qctool ### PHWE ≥ 1 × 10-7

# Step 4: GRM calculation (EUR & EAS)
for i in {1..50}
  do
  gcta  --mbfile genetic_file_list --keep qualified_sample --make-grm-part 50 $i --thread-num 50 --out genetic_grm 
done

cat genetic_grm.*.grm.id > genetic.grm.id
cat genetic_grm.*.grm.bin > genetic.grm.bin
cat genetic_grm.*.grm.N.bin > genetic.grm.N.bin

gcta --grm genetic --make-bK-sparse 0.05 --out genetic_sp_grm



# Step 5: GRM calculation (non-EUR)
Rscript genesis_make.R
gcta  --grm-gz grm_genesis --make-grm --out grm_gcta  --thread-num 20

gcta --grm grm_gcta --make-bK-sparse 0.05 --out grm_sp_grm

