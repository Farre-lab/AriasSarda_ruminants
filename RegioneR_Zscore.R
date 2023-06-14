library(ChIPpeakAnno)
library(rtracklayer)
library(ensembldb)
library(ChIPseeker)
library(GenomicFeatures)
library(regioneR)

# File containing the coordinates of the EBRs and HSBs 
# Structure: chromosome start end (only chromosome number: 1 (not chr1))
tads <- read.table("coordinates_EBRs.txt", sep="", stringsAsFactors = F, header=F)

# File containing the coordinates of: all cattle genes, housekeeping genes, and GC content specific 
# Structure: chromosome start end (only chromosome number: 1 (not chr1))
chip <- read.table("coordinates_cattle_genes.txt", sep="", stringsAsFactors = F, header=F)

# Cattle chromosome size
# Structure: chromosome start end (only chromosome number: 1 (not chr1))
genome <- toGRanges("Cattle_sizes.txt")

pt <- permTest(A=tads, B=chip, ntimes=1000, randomize.function=circularRandomizeRegions,count.once=TRUE,
               evaluate.function=numOverlaps, genome=genome, verbose = TRUE,  force.parallel=FALSE, per.chromosone=FALSE)

lz <- localZScore(pt=pt, A=tads, B=chip)


# plot RegioneR
plot(pt)

#plot local Z score
plot(lz)
