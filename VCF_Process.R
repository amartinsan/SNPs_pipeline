#Read files to then import them to R 
library(data.table)
library(vcfR)
library(SNPRelate)
#make a function to import a .vcf file 

importVCF <- function(vcffile, colClasses = "character", header = T, ...){
  fread(input = vcffile, colClasses = "character", header = T, ...)
}
#Function to ge random sting for later ID subsitution
#Can also be numbers 
#vals=seq(1:10000000)
RandomID <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#Read all the files to import them all in one list
#####It is important to setwd() on the carpet with all the files#
filenames <- list.files(getwd(), pattern="*.vcf", full.names=TRUE)
list_VCFfiles <- lapply(filenames, importVCF)

#Merge the list of vCF files
ALLsampleVCF=merged.data.frame = Reduce(function(...) merge(..., all=T), list_VCFfiles)

#Substitute the "." in ID to NA so we can change it, #####GENLIGHT CANNOT USE NOT UNIQUE IDs########
Original_all_samplesVCF = ALLsampleVCF
ALLsampleVCF$ID=sub(x = ALLsampleVCF$ID,".",replacement = NA,fixed = T)
vals=seq(1:10000000)

#is.na function its perfect to replace IDs with a random number
ALLsampleVCF$ID[is.na(ALLsampleVCF$ID)]= sample(x = vals,size = sum(is.na(ALLsampleVCF$ID)),replace = T)

#Use dplyr to filter duplicates

library(dplyr)
Dataframe_vCF <- ALLsampleVCF %>% arrange(ID) %>%filter(duplicated(ID) == FALSE)

#Both values have to be equal
length(Dataframe_vCF$ID)
length(unique(Dataframe_vCF$ID))

#Write the new data frame as a  VCF file. It needs to have metadta
#####Kinda tricky but you can add the metadata with a conection#####

con <- file("merge_genotype.vcf", open="wt")
writeLines(paste("##fileformat=VCFv4.1"), con)
write.table(Dataframe_vCF,con,quote=F,sep = "\t",row.names = F)
close(con)  
###Now We have the merged genotype A .vcf file with all samples#####
######Save data for storage and working after######
save(list_VCFfiles,Original_all_samplesVCF,Dataframe_vCF,importVCF, file="VCF.rda")


################
######################################ANALYSIS OF VCF samples#########################################
#################


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyr)
library(parallel)

Amphora.VCF <- read.vcfR("merge_genotype.vcf")
#Read the Superpopulation codes

coordination=read.table("C:/Users/AMS/Desktop/Amphora/DataEngineering-GenomicsChallenge-Jun2022/Challenge/coordination.txt", sep = "\t", header = TRUE)
#Convert vcf to genlight
gl.rubi <- vcfR2genlight(Amphora.VCF)

#we need to modify the pop.data for better use

gtnames=colnames(Amphora.VCF@gt)
df=data.frame(gtnames,row.names = NULL)
df=data.frame(df[-1,])
colnames(df) ="UUID"
row.names(df) = df$UUID
row.names(coordination) = coordination$UUID
population=merge(df,coordination,by =0,all = T)
population.data=population[,c(2,4)]
population.data$Superpopulation.code = population.data$Superpopulation.code %>% replace_na("NoCODE")
#One UUID is not in the samples
pop.data=population.data[-765,]
write.table(population.data,quote=F,sep = "\t",row.names = F,file = "Sample_UUIDNAMES.txt")

#Add the population data to the genlight file

pop(gl.rubi) <- pop.data$Superpopulation.code
#Define the plody
ploidy(gl.rubi) <- 2
#Save files
save(gl.rubi,pop.data,file="DATA_for_VCFr_analysis.rda")
#The vcfr file is big, save apart. Sometimes you only need the genlight file.
save(Amphora.VCF,file="VCFrfile.rda")



