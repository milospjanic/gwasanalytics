# bed2GwasCatalogBinomialGgplot

The P-values will be computed using binomial cumulative distribution function b(x;n,p) in R (dbinom function). We set the parameter n equal to the total number of GWAS SNPs in a particular GWAS phenotype. Parameter x was set to the number of GWAS SNPs for a given GWAS phenotype that overlap input regions and parameter p was set to the fraction of the uniquely mappable human hg19 genome (calculated with subscript) that is localized in the input regions. Calculated binomial p-value equals the probability of having x or more of the n test genomic regions in the input domains given that the probability of that occurring for a single genomic location is p. 

This script will connect to GWAS Catalog and download the entire data set, create GWAS Catalog bed file, and parse and uniq according to the N-1 input arguments. Argument provided to the bash script **runcategories.sh** should be a bed file that will be used to intersect parsed bed files from the GWAS Catalog. Script **runcategories.sh** contains the names of tested GWAS categories, modify it if you want different ones. Script **runcategories.sh** will then call **bed2GwasCatalogBinomialGgplot.sh**. Number of overlaps is reported and initial number of entries in parsed files. Finally, **bed2GwasCatalogBinomialGgplot.sh** will create an R script that will be executed to calculate binomial p-values for each overlap. Intermediary files will be removed, except: GwasCatalog.bed (entire catalog in a bed file), \*gwascatalog.bed (parsed original files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.chrXY (parsed and uniqed files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.overlap (overlap with parsed files), GwasCatalog2Bed.sh (sciprt to download GWAS Catalog and convert to bed). Note that names of phenotypes in GWAS Catalog start with capital letter but then next word is with small letter. **That is why we enabled case insensitive search in this script.**

#Usage
<pre>
chmod 775 ./bed2GwasCatalogBinomial.sh
./bed2GwasCatalogBinomial "Coronary artery" "Coronary heart" "Bipolar disorder" "Feminism" L2_TCCTGAGC_L002_peaks.bed 
</pre>

#Dependencies 
Rscript, bedtools (needs to be in $PATH)

#Output
<pre>
