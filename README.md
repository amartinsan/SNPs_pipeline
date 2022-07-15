 # Deliverable for Amphora Health Challenge
 Code used for processing data from DataEngineering-GenomicsChallenge-Jun2022.zip.
 
 [Challenge](https://github.com/AdrianMtz-Santana/AmphoraHealth_Bioinformatic_Challenge2022_deliverable/blob/main/Coding_Challenge%20Instructions_Data%20Engineer_Genomics_Aug2022.pdf)
 
 Samples were processed using some bash functions to change the 382 raw files in tab-separated values to VCF.
 
 - Using the script: [csv_to_vcf.sh](https://github.com/AdrianMtz-Santana/AmphoraHealth_Bioinformatic_Challenge2022_deliverable/blob/397eaf7bb3acfd12c30d0cd7b15e53ce2e3eb538/csv_to_vcf.sh)
 
 code goes like this:
   
   ### First, we unzip the data and make some symbolic links to our working directory

    mkdir AmphoraChallenge
    cd AmphoraChallenge
    ln -s ~/DataEngineering-GenomicsChallenge-Jun2022/Challenge/Challenge\ Samples/Challenge\ Samples/* .
   ### The .csv files lack some processing to transform to a .vcf file.	
   ### First, the structure has to be arranged for better manipulation.
   ### This can be done in R or Python, but with simple BASH and AWK without importing to Python or R before Clustering or VCF analysis.
   ### Files have some columns: CHROM, POS, REF, ALT, and part of FORMAT in GT.
   
    for file in *.csv ;
    do 	
    #sed has a confusing  syntax; more lines of sed instead of a long one with all the changes gets clearer 
    
	  sed -i 's/,REF/#CHROM;POS,REF/' $file ; 
	  sed -i 's/,ALT,ALT/,ALT/' $file ;
	  sed -i 's/;/\t/' $file ;
	  sed -i 's/,/\t/' $file ;
	  sed -i 's/,/|/4g' $file ;
	  sed -i 's/"//' $file ;
	  sed -i 's/"//' $file ;
	  sed -i 's/,/\t/' $file ;
	  sed -i 's/,/\t/2g' $file ;
	  sed -i 's/ALT,/ALT\t/' $file ;
   
    #We are missing ID, QUAL, FILTER, INFO and FORMAT, columns we can add them with AWK
    
       awk 'BEGIN{ FS=OFS="\t" } {$2 = $2 FS (NR==1? "ID" : ".") }1' $file > tmp && mv tmp $file ;
	   awk 'BEGIN{ FS=OFS="\t" } {$5 = $5 FS (NR==1? "QUAL": ".") }1' $file > tmp && mv tmp $file ;
	   awk 'BEGIN{ FS=OFS="\t" } {$6 = $6 FS (NR==1? "FILTER" : "PASS") }1' $file > tmp && mv tmp $file ;
	   awk 'BEGIN{ FS=OFS="\t" } {$7 = $7 FS (NR==1? "INFO" : "AC=1;AN=2") }1' $file > tmp && mv tmp $file ;
	   awk 'BEGIN{ FS=OFS="\t" } {$8 = $8 FS (NR==1? "FORMAT" : "GT") }1' $file > tmp && mv tmp $file ;
    
    #Head just to check structure in logout 
	  
    head $file 
    
    #The following is to add metadata to maintain format consistency

	  sed -i 1i"##" $file 
	  sed -i 1i"##File was a .csv raw data, transformed to a .vcf file" $file 
	  sed -i 1i"##fileformat=VCFv4.1" $file ;
	  
    done 1> CSVtoVCF_logout.txt 2> CSVtoVCF_error.txt
   
### Change file extension 
	for file in *.csv; 
	do 
	mv -- "$f" "${f%.csv}.vcf"
	done

    
 ## The rest of the analysis was made using R 
 - The code is in: [VCF_Process.R](https://github.com/AdrianMtz-Santana/AmphoraHealth_Bioinformatic_Challenge2022_deliverable/blob/397eaf7bb3acfd12c30d0cd7b15e53ce2e3eb538/VCF_Process.R)
 
 To process the data i had to read A LOT of online guides 🥴 and watch many many videos.
 
 most of the help came from: 
 - https://grunwaldlab.github.io/Population_Genetics_in_R/index.html
 - https://github.com/cmcouto-silva?tab=repositories 
 - https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#overview
 
 ## List of Libraries to install 
 - data.table
 - vcfr
 - SNPRelate 
 - dplyr
 - poppr
 - ape 
 - RColorBrewer
 - tidyr
 - parllell
 - SeqArray
 - ggploy 
 - 
 
 
 
 ## Results for the clustering/ordination of all provided samples are shown in Results folder
 
 -Also the PNG graphs and merged_genotype.vcf file with all samples.
 
 
 ## For the file to take up new data, check the folder Feed_newDATA
 
 The code to work requires to run in the same folder as the file intented to parse to the large merge_genotype.vcf 
 also with the PCA-vcfa-data.rda file. 
 
 
 ## Perspectives
 
- Other ordination/clustering techniques UMAP or t-SNE techniches are kinda new techniques that can cluster samples in an efficient way. 

Althought that does not mean that PCA, PCoA or NMDS classical ordinations cant work well for some data. 

For example when you are testing stuff like dummy data

- This can be scalable in a server with Nextflow, Docker or API (Flask)

- Betwern VCFr and SeArray/SNPrelate libraries a testing to compare wich is faster or efficient for big ammounts of data and samples should be done.



