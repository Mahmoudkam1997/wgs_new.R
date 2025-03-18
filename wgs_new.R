rm(list=ls())
library(ShortRead)       # For handling FASTQ files
library(Biostrings)      # For sequence operations
library(Rsubread)        # For alignment
library(VariantAnnotation) # For variant calling
library(GenomicAlignments) # For BAM file manipulation
library(tidyverse)
library(AnnotationHub)
fastq_wgs <- "E:/fastq WGS/"
reads <- readFastq(fastq_wgs)
############################################################
qu <- quality(reads)
summary(qu)
####################################################
qa <- qa(fastq_wgs, type = "fastq")
report(qa, dest="quality_report")
#########################################################
trimmed_reads <- trimTails(reads, 5, "1", consecutive = TRUE)
filtered_reads <- reads[width(reads) > 50]
####################################################################
writeFastq(filtered_reads, "filtered_reads.fastq.gz", compress=TRUE)
#####################################################################
reference_genome <- "E:/sequence.fasta"
buildindex(basename = "ref_index", reference = reference_genome)
######################################################################
# Align reads
align(index = "ref_index",
      readfile1 = "filtered_reads.fastq.gz",
      output_file = "aligned_reads.bam",
      nthreads = 4)
##############################################################
library(Rsamtools)

# Sort BAM file
sorted_bam <- "E:/aligned_reads.sorted.bam.bam"
sortBam("E:/aligned_reads.bam", sorted_bam)

# Index BAM file
indexBam(sorted_bam)
################################################################
#variant calling in Galaxy 
############################################################
#variant annotation
vcf <- readVcf("E:/Galaxy15-[bcftools mpileup on data 12 and data 14].vcf", 
               genome = "E:/sequence.fasta")
info(vcf)
#################################################################
#Annotate variants another way 
# Load the genome annotation
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("E:/Escherichia_coli_gca_003856675.ASM385667v1.60.gff3")

# Extract genes
genes <- genes(txdb)

# Find overlaps between variants and genes
overlaps <- findOverlaps(rowRanges(vcf), genes)

# View annotated variants
annotated_variants <- genes[subjectHits(overlaps)]
annotated_variants
##############################################################
BiocManager::install("Gviz")
library(Gviz)
# Load BAM file
bamTrack <- AlignmentsTrack("E:/aligned_reads.sorted.bam.bam", 
                            genome = "txdb")

# Create a genome axis track
genomeAxis <- GenomeAxisTrack()

# Plot the tracks
plotTracks(list(genomeAxis, bamTrack), from = 100000, to = 101000)
#####################################################################
library(clusterProfiler)
BiocManager::install("org.EcK12.eg.db")  # For *E. coli*
# Convert gene IDs to KEGG
library(org.EcK12.eg.db)
kegg_genes <- bitr(annotated_variants$gene_id, fromType = "GENENAME",
                   toType = "PATH", OrgDb = org.EcK12.eg.db)

entrez_ids <- bitr(annotated_variants$gene_id, fromType = "GENENAME",
                   toType = "ENTREZID", OrgDb = "org.EcK12.eg.db")
# Perform KEGG enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                          organism = "eco")
head(kegg_enrich)
dotplot(kegg_enrich, showCategory = 20)
emapplot(kegg_enrich)
cnetplot(kegg_enrich, showCategory = 10, circular = TRUE, colorEdge = TRUE)



