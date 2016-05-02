# gwasanalytics

Gwas Analytics consists of a series of scripts that will analyze Gwas Catalog phenotypes and the input bed file from the ChIP-Seq experiment (or any other experiment that involves patritioning of the genome into compartments, e.g. ATAC-Seq, DNase-Seq etc).

#bed2GwasCatalogBinomialGgplot
This script will calculate binomial p-values and fold change of genomic overlaps between the input bed file and various GWAS categories.

#bed2GwasCatalogBinomialMod1Ggplot
The modification of the binomial test includes calculation of the binomial p-value using modified probability of the base nucleotide to localize in the tested region. In the standard binomial test, one uses the probability of the base nucleotide to be localized in the tested input regions of the bed file, which corresponds to the percent of input bed in the human genome. This approach is robust and less sensitive towards the larger (but not smaller) sizes of the tested GWAS category, however it does not resolve well the p-values of various categories, which you can see in the examples provides. In addition, calculated P-value will approximately correspond to the fold change (number of overlaps/total number of GWAS variants), due to the fact in addition to number of overlaps/total number of GWAS variants (that correspond to the fold change) the only parameter used to calculate binomial p-value is the probability of a base to be in the input region and which is constant for each tested GWAS category. Hence, this is the reason the p-value/fold change plot will look almost diagonal in every example.


Therefore, we modified the calculation of the probability to include the probability of base nucleotide to be present in the input regions that correspond to the tested GWAS category. This modification will better resolve p-values and will be less sensitive to the inflated p-values of smaller datasets, however the modification itself inflates p-values, which is something to consider when interpreting the results. One example of the better resolution of p-values is the example of open chromatin regions in ENCODE Glioblastoma cell line. When tested with the standard binomial test the brain related categories such as Schizophrenia and Bipolar disorder do not emerge tot the top of the list. Surprisingly, the top categories are not brain related and are relatively small datasets (CardiogramplusC4D, Diastolic blood pressure, Coronary heart disease). When modified binomial is applied, we obtained better resolution of GWAS categories, with the brain related categories, Schizophrenia and Bipolar disorder, being the top ones.

**Example of standard binomial test - Glioblastoma ENCODE cell line**. GWAS categories did not resolve well according to the p-values. Negative log p-value and fold change are proportional and the graph is nearly a diagonal.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/Glioblastoma.ENCODE.binomial.png)

**Example of modified binomial test - Glioblastoma ENCODE cell line**. We can see better resolution of GWAS categories by their p-values. Note that modified binomial gives brain related categories (Schizophrenia and Bipolar disorder) as top ones, as expected for a cell line with the brain origin, such as Glioblastoma ENCODE cell line.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/Glioblastoma.ENCODE.binomial.mod.png)

Another example of the better separation of p-values is the example of ATAC-Seq open chromatin regions in CD4+ T- cells. When tested with the standard binomial test the immune-related categories do not emerge to the top of the list, while the top categories are not immune- related and are again either relatively small or unrelated datasets (CardiogramplusC4D, Coronary Heart disease). When modified binomial is applied, we obtained better resolution of GWAS categories, with the immune- category, Rheumatoid arthritis, being the top category and another immune- category, Crohns disease, being among the top categories. Especially, this is relevant as it is known that Rheumatoid arthritis is a chronic inflammatory syndrome with a central role of CD4+ T-cells and is characterized with aberrant pathways of CD4+ T-cell activation.

**Example of standard binomial test - CD4+ T-Cell ATAC-Seq GEO:GSE60682** GWAS categories did not resolve well according to the p-values. Negative log p-value and fold change are proportional and the graph is nearly a diagonal.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/CD4TCell.GSE60682.binomial.png)

**Example of modified binomial test - CD4+ T-Cell ATAC-Seq GEO:GSE60682** We can see better resolution of GWAS categories by their p-values. Note that modified binomial gives immune- related categories (Rheumatoid arthritis, Crohns disease) as top ones, as expected for a CD4 T cell, one of the main components of the immune system.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/CD4TCell.GSE60682.binomial.mod.png)


Here is an example of the better separation of p-values in the case of ATAC-Seq open chromatin regions in GM12878 cell line, a lymphoblastoid cell line produced from the blood (B-lyphocytes) of a female donor with northern and western European ancestry. When tested with the standard binomial test the top category is relatively small and unrelated Diastolic blood pressure dataset. When modified binomial is applied, we obtained better resolution of GWAS categories, with the immune- category, Crohns disease, being the top category and other immune- categories, Multiple sclerosis, Asthma and Ulcerative colitis, being among the top categories, which is consistent with the fact that this cell line orginates from an immune cell. In addition we obtained much higher ranking for some types of cancers, e.g. breast cancer, which is in line with the fact that these cells originate from a (lymphoblastoid) cancer.

**Example of standard binomial test - GM12878 cell line ATAC-Seq GEO:GSE47753** GWAS categories did not resolve well according to the p-values. Negative log p-value and fold change are proportional and the graph is nearly a diagonal.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/GM12878.GSE47753.binomial.png)

**Example of modified binomial test - GM12878 cell line ATAC-Seq GEO:GSE47753** We can see better resolution of GWAS categories by their p-values. Note that modified binomial gives immune- related categories (Crohns disease, Multiple sclerosis, Asthma and Ulcerative colitis) as top ones, which is consistent with the immune- origin of GM12878.

![ScreenShot](https://github.com/milospjanic/IntegrativeFunctionalGenomics/blob/master/GM12878.GSE47753.binomial.mod.png)
