# The-genetic-architecture-of-brainstem-structures
Code for research paper: The genetic architecture of brainstem structures

# Imaging preprocess
Please see 1_imaging_prep.sh
Software
  1. Brainstem substructure segmentation: FreeSurfer v7.0.0 (https://surfer.nmr.mgh.harvard.edu/)
  2. Harmonzationt: Combat Harmonization (https://github.com/Jfortin1/ComBatHarmonization)
  3. Quantile normalization: PLINK (https://www.cog-genomics.org/plink/1.9/)
# Genetic preprocess
Please see 2_genetic_prep.sh
Software
  1. Ancestry: SNPWeights (https://hsph.harvard.edu/research/price-lab/software/)
  2. Genetic PC calculatation: ABCD_GRM.r (https://github.com/robloughnan/ABCD_GeneticPCs_and_Relatedness)
  3. Marker-level QC
    a. PLINK (https://www.cog-genomics.org/plink/1.9/)
    b. PLINK2 (https://www.cog-genomics.org/plink/2.0/)
    c. QCtool (https://www.chg.ox.ac.uk/~gav/qctool_v2/)
  5. GRM calculation
    a. EUR & EAS: GCTA-GRM (https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM)
    b. non-EUR: GENESIS (https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)
# The genetic architecture of brainstem structures
Please see 3_main.sh
Software
  1. GWAS: GCTA-fastGWA (https://yanglab.westlake.edu.cn/software/gcta/#fastGWA)
  2. GWAS meat analyses: METAL (https://genome.sph.umich.edu/wiki/METAL_Documentation)
  3. Pooling: PLINK (https://www.cog-genomics.org/plink/1.9/)
  4. h2
       a. EUR & EAS: LDSC (https://github.com/bulik/ldsc)
       b. EUR+EAS: cov-LDSC (https://github.com/immunogenomics/cov-ldsc)
  5. cross-ancestry genetic correlation: Popcorn (https://github.com/brielin/Popcorn)
 # Biological mechanisms
 Please see 4_post_gwas.sh
 Software
   1. Fine-mapping: PAINTOR (https://github.com/gkichaev/PAINTOR_V3.0)
   2. Annotation: FUMA (https://fuma.ctglab.nl/)
   3. Gene priotization: FLAMES (https://github.com/Marijn-Schipper/FLAMES)
        a. MAGMA (https://cncr.nl/research/magma/)
        b. PoPs (https://github.com/FinucaneLab/pops)
   5. Pathway enrichment analyses: g:Profiler (https://biit.cs.ut.ee/gprofiler/gost)
# Overlap with other phenotypes
Please see 5_pheno.sh
Software
  1. PheWAS: FinnGen-MVP-UKBB (https://mvp-ukbb.finngen.fi/)
  2. genetic correlation: LDSC (https://github.com/bulik/ldsc)
  3. genetic colocalization: coloc (https://chr1swallace.github.io/coloc/)
  4. conjFDR: pleioFDR (https://github.com/precimed/pleiofdr)
