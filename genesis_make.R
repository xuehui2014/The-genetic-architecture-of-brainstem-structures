library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(optparse)
snpgdsBED2GDS(bed.fn = "genotype.bed", bim.fn = "genotype.bim", fam.fn = "genotype.fam", out.gdsfn = "mygenotype.gds")
geno <- GdsGenotypeReader(filename = "mygenotype.gds")
myenoData <- GenotypeData(geno)
gds <- snpgdsOpen(paste(genotype,".gds",sep=""))
# Kinship Matrix 
ibd.robust = snpgdsIBDKING(gds)
ibd_king <- snpgdsIBDKING(gdsfile, type="KING-robust", verbose=TRUE)
row.names(KINGmat) = ibd.robust$sample.id
colnames(KINGmat) = ibd.robust$sample.id
# PCA-AiR
ABCD_geno <- GdsGenotypeReader(filename = gds.fn)
# create a GenotypeData class object
ABCD_genoData <- GenotypeData(ABCD_geno)
ABCD_genoData

print('Generating pcair')
mypcair <- pcair(ABCD_genoData, kinobj = KINGmat, divobj = KINGmat, 
                 snp.include = pruned)

pcair_file = paste0(fn_prefix, '_pcair.RData')
save(mypcair, file=pcair_file)
print(paste0('Saved pcair ', pcair_file))

# Save unrelated individuals
unrel_file = paste0(fn_prefix, '_unrelateds.txt')
write.table(mypcair$unrels, unrel_file, quote=FALSE, col.names=FALSE, row.names=FALSE)

pcs = mypcair$vectors
colnames(pcs) = paste0('C', seq(dim(pcs)[2]))
pcair_file = paste0(fn_prefix, '_pcair.tsv')
write.table(pcs, pcair_file, quote=FALSE, sep='\t', col.names=NA)

# Generate loadings
loadings = pca_loadings(ABCD_geno@handler, kinobj = KINGmat, divobj = KINGmat, 
                 snp.include = pruned)
write.table(loadings, paste0(fn_prefix, '_pcair_loadings.tsv'), quote=FALSE, sep='\t', col.names=NA)

# PCA-Relate
pcrelate_file = paste0(fn_prefix, '_pcrelate.RData')

    print('Generating PC-Relate')
    ABCD_genoData <- GenotypeBlockIterator(ABCD_genoData, snpInclude=pruned)
    mypcrelate <- pcrelate(ABCD_genoData, pcs = mypcair$vectors[,1:2], 
                        training.set = mypcair$unrels)
    save(mypcrelate, file=pcrelate_file)
    print(paste0('Saved PC-Relate ', pcrelate_file))
# Create GRM 
print('Generating GRM')
grm = pcrelateToMatrix(mypcrelate)
grm_mat = format(as.matrix(grm), digits=8)
grm_out = paste0(fn_prefix, '_GRM.tsv')
write.table(grm_mat, grm_out, quote=FALSE, sep='\t', col.names=NA)
print(paste0('Saved GRM to ', grm_out))




