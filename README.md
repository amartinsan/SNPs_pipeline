 # Deliverabe for Amphora Health Challenge
 Code used for processing data from DataEngineering-GenomicsChallenge-Jun2022.zip.
 
 Samples where proccesed by using some bash functions to change the 382  raw files in tab-separated values to VCF.
 - Using the script: csv_to_vcf.sh
 
 code goes like this:
   
   ### First we unzip the data and make some simbolic links to our working directory

    mkdir AmphoraChallenge
    cd AmphoraChallenge
    ln -s ~/DataEngineering-GenomicsChallenge-Jun2022/Challenge/Challenge\ Samples/Challenge\ Samples/* .
   ### The .csv files lack some proccesing to transform to a .vcf file.	
   ### First the structure has to be arranged to better manipulation.
   ### This can be done in R or Python, but also with simple BASH and AWK without the need of importing to Python or R before Clustering or VCF analysis.
   ### Files already have some columns: CHROM,POS,REF,ALT, and part of FORMAT in GT.
   
    for file in *.csv ;
    do 	
    #sed has a confunsing syntax, more lines of sed instead of a long one with all the changes gets clearer 
    
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
