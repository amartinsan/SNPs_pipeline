library(data.table)
library(vcfR)
library(SNPRelate)
library(SeqArray)
library(dplyr)

load("PCA-vcfa-data.rda")
importVCF <- function(vcffile, colClasses = "character", header = T, ...){
  fread(input = vcffile, colClasses = "character", header = T, ...)
}
RandomID <- function(n = 5000) {
  a <- do.call(paste0, replicate(sample(1:6, 1), sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
userfile <- importVCF(fn)
#Merge all samples to new file
forjoin=colnames(userfile)
forjoin=forjoin[1:9]
ALLsampleVCF=merge(userfile,Dataframe_vCF,by = forjoin, all = T)
#Replace "." for "Na" ind  ID column
ALLsampleVCF$ID=sub(x = ALLsampleVCF$ID,".",replacement = NA,fixed = T)
set.seed(sample(1:15, 1))
vals=RandomID(10030)
#Filter duplicates
#is.na function its perfect to replace IDs with a rando string
ALLsampleVCF$ID[is.na(ALLsampleVCF$ID)]= sample(x = vals,size = sum(is.na(ALLsampleVCF$ID)),replace = T)
#Use dplyr to filter duplicates
#merged_VCF_DF <- ALLsampleVCF %>% arrange(ID) %>%filter(duplicated(ID) == FALSE)
merged_VCF_DF=Dataframe_vCF2=unique(setDT(ALLsampleVCF)[order(ID)], by = "ID")


#####Save table with new sample#
con <- file("NEWSAMPLE_merge_genotype.vcf", open="wt")
writeLines(paste("##fileformat=VCFv4.1"), con)
write.table(merged_VCF_DF,con,quote=F,sep = "\t",row.names = F)
close(con) 


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
#IMPORT data
vcf.fn <- "NEWSAMPLE_merge_genotype.vcf"
snpgdsVCF2GDS(vcf.fn, "GDS_file.gds", method="biallelic.only")
seqVCF2GDS(vcf.fn, "GDS_file.gds")
genofile <- seqOpen("GDS_file.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))
#Calculate
pca <- snpgdsPCA(genofile, snp.id=snpset.id, autosome.only = F,remove.monosnp = F,bayesian = T)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#add a row to superpopulation data
popsample[nrow(popsample) + 1,] = c(colnames(Testfile[,10]),colnames(Testfile[,10]),"EXTRASAMPLE")
#get graph
tablePCA <- data.frame(sample.id = pca$sample.id,pop = factor(popsample$Superpopulation.code)[match(pca$sample.id, sample.id)],PC1 = pca$eigenvect[,1],PC2 = pca$eigenvect[,2],stringsAsFactors = FALSE)

#Colors vector "paired can get only 12 
cols <- brewer.pal(n =nlevels(as.factor(popsample$Superpopulation.code)), name = "Paired")
plotPCA=ggplot(tablePCA, aes(x=PC1, y=PC2, colour=pop)) + geom_point(size=2) + scale_color_manual(values = cols) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+theme_bw() + ggtitle("PCA with all samples") + theme(legend.position = "bottom", legend.title = element_blank(), axis.title = element_text(size = 17), axis.text = element_text(size = 14), legend.text = element_text(size = 15))
#Export plot
jpeg(file="PCA.jpeg")
ggplot(tablePCA, aes(x=PC1, y=PC2, colour=pop)) + geom_point(size=2) + scale_color_manual(values = cols) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+theme_bw() + ggtitle("PCA with all samples") + theme(legend.position = "bottom", legend.title = element_blank(), axis.title = element_text(size = 17), axis.text = element_text(size = 14), legend.text = element_text(size = 15))
dev.off()
#Identity-By-State Analysis
ibs <- snpgdsIBS(genofile)
ibsmax=data.frame(ibs$ibs)
#Substitute NaN to 0
ibsmax[is.na(ibsmax)] <- 0
ibs$ibs = as.matrix(ibsmax)
loc <- cmdscale(1 - ibs$ibs, k = 2)
loz=as.data.frame(loc)
race <- as.factor(popsample$Superpopulation.code)
MDSPLOT=ggplot(loz, aes(x=V1, y=V2, colour=race)) + geom_point(size=2) +scale_color_manual(values=cols)+geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() +theme(legend.title = element_blank(), axis.title = element_blank(), axis.text = element_text(size = 14), legend.text = element_text(size = 15))
jpeg(file="PCA.jpeg")
dev.off()
#save plots data
save(MDSPLOT,plotPCA,file = "DATA_plots.rda")
