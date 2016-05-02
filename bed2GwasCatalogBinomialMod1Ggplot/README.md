#bed2GwasCatalogBinomialMod1Ggplot

This script is a modification of the bed2GwasCatalogBinomial and calculates binomial p-value for genomics overlaps using the following criteria. The P-values were computed using binomial cumulative distribution function b(x;n,p) in R (dbinom function). We set the parameter n equal to the total number of GWAS SNPs in a particular GWAS phenotype. Parameter x was set to the number of GWAS SNPs for a given GWAS phenotype that overlap input regions and parameter p was set to the fraction of the uniquely mappable human hg19 genome (calculated with subscript) that is localized in the input regions and contains assessed GWAS phenotype SNPs. Calculated binomial p-value equals the probability of having x or more of the n test genomic regions in the open chromatin domain given that the probability of that occurring for a single GWAS genomic location is p. Note the p-values may be inflated depending on your input files. 

This script will connect to GWAS Catalog and download the entire data set, create bed file, and parse and uniq according to the N-1 input arguments. Last argument provided to the bash script should be a bed file that will be used to intersect parsed bed files from the GWAS Catalog. Number of overlaps is reported and initial number of entries in parsed files.  Finally, it will create an R script that will be executed to calculate binomial p-values for each overlap. Intermediary files will be removed, except: GwasCatalog.bed (entire catalog in a bed file), \*gwascatalog.bed (parsed original files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.chrXY (parsed and uniqed files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.overlap (overlap with parsed files, GWAS SNP positions), \*gwascatalog.bed.cut.sort.uniq.overlap.input.int.cut (overlap with parsed files, input regions), GwasCatalog2Bed.sh (sciprt to download GWAS Catalog and convert to bed). Note that names of phenotypes in GWAS Catalog start with capital letter but then next word is with small letter. **That is why we enabled case insensitive search in this script.**

#Usage
<pre>
chmod 775 runcategories.sh 
./runcategories.sh GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed 
</pre>

#Dependencies 
Rscript, bedtools (needs to be in $PATH)

#Output
<pre>
chmod 775 runcategories.sh 
./runcategories.sh GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed

GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed
--2016-05-02 00:52:22--  http://www.genome.gov/admin/gwascatalog.txt
Resolving www.genome.gov (www.genome.gov)... 156.40.242.24
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:80... connected.
HTTP request sent, awaiting response... 301 Moved Permanently
Location: https://www.genome.gov/admin/gwascatalog.txt [following]
--2016-05-02 00:52:22--  https://www.genome.gov/admin/gwascatalog.txt
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 10407265 (9.9M) [text/plain]
Saving to: `gwascatalog.txt'

100%[======================================================================================================================================================>] 10,407,265  1.30M/s   in 8.0s    

2016-05-02 00:52:31 (1.23 MB/s) - `gwascatalog.txt' saved [10407265/10407265]

Received: Bipolar disorder
Done: Bipolar disorder
Received: Schizophrenia
Done: Schizophrenia
Received: Crohn's disease
Done: Crohn's disease
Received: Rheumatoid arthritis
Done: Rheumatoid arthritis
Received: Bone mineral density
Done: Bone mineral density
Received: Ulcerative colitis
Done: Ulcerative colitis
Received: Coronary Heart
Done: Coronary Heart
Received: Type 2 diabetes
Done: Type 2 diabetes
Received: Alzheimer's disease
Done: Alzheimer's disease
Received: Body mass index
Done: Body mass index
Received: Multiple sclerosis
Done: Multiple sclerosis
Received: Breast cancer
Done: Breast cancer
Received: Asthma
Done: Asthma
Received: Triglycerides
Done: Triglycerides
Received: Parkinson's disease
Done: Parkinson's disease
Received: LDL cholesterol
Done: LDL cholesterol
Received: Blood pressure
Done: Blood pressure
Received: Coronary Artery
Done: Coronary Artery
Received: Coronary artery calcification
Done: Coronary artery calcification
Received: Prostate cancer
Done: Prostate cancer
Received: Liver enzyme levels
Done: Liver enzyme levels
Received: Diastolic blood pressure
Done: Diastolic blood pressure
Received: Hypertension
Done: Hypertension
Received: HDL cholesterol
Done: HDL cholesterol
Received: QRS duration
Done: QRS duration
Received: Primary biliary cirrhosis
Done: Primary biliary cirrhosis
Received: BMI
Done: BMI
Received: Type 1 diabetes
Done: Type 1 diabetes
Received: Longevity
Done: Longevity
Received: Waist-hip ratio
Done: Waist-hip ratio
Received: Inflammatory biomarkers
Done: Inflammatory biomarkers
Received: Lipid metabolism phenotypes
Done: Lipid metabolism phenotypes
Received: Kawasaki disease
Done: Kawasaki disease
Received: Pancreatic cancer
Done: Pancreatic cancer
Received: C-reactive protein
Done: C-reactive protein
Received: Celiac disease
Done: Celiac disease
Received: Lupus
Done: Lupus
Received: Myocardial infarction
Done: Myocardial infarction
Received: Insulin resistance
Done: Insulin resistance
Received: FEV1
Done: FEV1
Received: Fasting plasma glucose
Done: Fasting plasma glucose
Received: Menarche
Done: Menarche
Received: CardiogramplusC4D
Done: CardiogramplusC4D
Gwas Catalog number of SNP-phenotype associations:
18951 GwasCatalog.bed
Gwas Catalog number of SNP-phenotype associations per category:
Phenotype: Bipolar disorder
444 Bipolar disorder.gwascatalog.bed
Phenotype: Schizophrenia
413 Schizophrenia.gwascatalog.bed
Phenotype: Crohn's disease
228 Crohn's disease.gwascatalog.bed
Phenotype: Rheumatoid arthritis
369 Rheumatoid arthritis.gwascatalog.bed
Phenotype: Bone mineral density
245 Bone mineral density.gwascatalog.bed
Phenotype: Ulcerative colitis
152 Ulcerative colitis.gwascatalog.bed
Phenotype: Coronary Heart
143 Coronary Heart.gwascatalog.bed
Phenotype: Type 2 diabetes
318 Type 2 diabetes.gwascatalog.bed
Phenotype: Alzheimer's disease
226 Alzheimer's disease.gwascatalog.bed
Phenotype: Body mass index
265 Body mass index.gwascatalog.bed
Phenotype: Multiple sclerosis
215 Multiple sclerosis.gwascatalog.bed
Phenotype: Breast cancer
253 Breast cancer.gwascatalog.bed
Phenotype: Asthma
166 Asthma.gwascatalog.bed
Phenotype: Triglycerides
180 Triglycerides.gwascatalog.bed
Phenotype: Parkinson's disease
105 Parkinson's disease.gwascatalog.bed
Phenotype: LDL cholesterol
166 LDL cholesterol.gwascatalog.bed
Phenotype: Blood pressure
266 Blood pressure.gwascatalog.bed
Phenotype: Coronary Artery
110 Coronary Artery.gwascatalog.bed
Phenotype: Coronary artery calcification
55 Coronary artery calcification.gwascatalog.bed
Phenotype: Prostate cancer
225 Prostate cancer.gwascatalog.bed
Phenotype: Liver enzyme levels
74 Liver enzyme levels.gwascatalog.bed
Phenotype: Diastolic blood pressure
36 Diastolic blood pressure.gwascatalog.bed
Phenotype: Hypertension
30 Hypertension.gwascatalog.bed
Phenotype: HDL cholesterol
213 HDL cholesterol.gwascatalog.bed
Phenotype: QRS duration
59 QRS duration.gwascatalog.bed
Phenotype: Primary biliary cirrhosis
38 Primary biliary cirrhosis.gwascatalog.bed
Phenotype: BMI
96 BMI.gwascatalog.bed
Phenotype: Type 1 diabetes
133 Type 1 diabetes.gwascatalog.bed
Phenotype: Longevity
50 Longevity.gwascatalog.bed
Phenotype: Waist-hip ratio
14 Waist-hip ratio.gwascatalog.bed
Phenotype: Inflammatory biomarkers
29 Inflammatory biomarkers.gwascatalog.bed
Phenotype: Lipid metabolism phenotypes
91 Lipid metabolism phenotypes.gwascatalog.bed
Phenotype: Kawasaki disease
20 Kawasaki disease.gwascatalog.bed
Phenotype: Pancreatic cancer
44 Pancreatic cancer.gwascatalog.bed
Phenotype: C-reactive protein
72 C-reactive protein.gwascatalog.bed
Phenotype: Celiac disease
90 Celiac disease.gwascatalog.bed
Phenotype: Lupus
174 Lupus.gwascatalog.bed
Phenotype: Myocardial infarction
13 Myocardial infarction.gwascatalog.bed
Phenotype: Insulin resistance
29 Insulin resistance.gwascatalog.bed
Phenotype: FEV1
9 FEV1.gwascatalog.bed
Phenotype: Fasting plasma glucose
17 Fasting plasma glucose.gwascatalog.bed
Phenotype: Menarche
82 Menarche.gwascatalog.bed
Phenotype: CardiogramplusC4D
52 CardiogramplusC4D.gwascatalog.bed
Gwas Catalog number of SNP-phenotype associations per category AFTER REMOVING DUPLICATES:
Phenotype: Bipolar disorder
411 Bipolar disorder.gwascatalog.bed.cut.sort.uniq
Phenotype: Schizophrenia
376 Schizophrenia.gwascatalog.bed.cut.sort.uniq
Phenotype: Crohn's disease
189 Crohn's disease.gwascatalog.bed.cut.sort.uniq
Phenotype: Rheumatoid arthritis
246 Rheumatoid arthritis.gwascatalog.bed.cut.sort.uniq
Phenotype: Bone mineral density
167 Bone mineral density.gwascatalog.bed.cut.sort.uniq
Phenotype: Ulcerative colitis
120 Ulcerative colitis.gwascatalog.bed.cut.sort.uniq
Phenotype: Coronary Heart
130 Coronary Heart.gwascatalog.bed.cut.sort.uniq
Phenotype: Type 2 diabetes
212 Type 2 diabetes.gwascatalog.bed.cut.sort.uniq
Phenotype: Alzheimer's disease
202 Alzheimer's disease.gwascatalog.bed.cut.sort.uniq
Phenotype: Body mass index
234 Body mass index.gwascatalog.bed.cut.sort.uniq
Phenotype: Multiple sclerosis
198 Multiple sclerosis.gwascatalog.bed.cut.sort.uniq
Phenotype: Breast cancer
196 Breast cancer.gwascatalog.bed.cut.sort.uniq
Phenotype: Asthma
161 Asthma.gwascatalog.bed.cut.sort.uniq
Phenotype: Triglycerides
117 Triglycerides.gwascatalog.bed.cut.sort.uniq
Phenotype: Parkinson's disease
87 Parkinson's disease.gwascatalog.bed.cut.sort.uniq
Phenotype: LDL cholesterol
101 LDL cholesterol.gwascatalog.bed.cut.sort.uniq
Phenotype: Blood pressure
171 Blood pressure.gwascatalog.bed.cut.sort.uniq
Phenotype: Coronary Artery
89 Coronary Artery.gwascatalog.bed.cut.sort.uniq
Phenotype: Coronary artery calcification
55 Coronary artery calcification.gwascatalog.bed.cut.sort.uniq
Phenotype: Prostate cancer
186 Prostate cancer.gwascatalog.bed.cut.sort.uniq
Phenotype: Liver enzyme levels
70 Liver enzyme levels.gwascatalog.bed.cut.sort.uniq
Phenotype: Diastolic blood pressure
34 Diastolic blood pressure.gwascatalog.bed.cut.sort.uniq
Phenotype: Hypertension
30 Hypertension.gwascatalog.bed.cut.sort.uniq
Phenotype: HDL cholesterol
141 HDL cholesterol.gwascatalog.bed.cut.sort.uniq
Phenotype: QRS duration
59 QRS duration.gwascatalog.bed.cut.sort.uniq
Phenotype: Primary biliary cirrhosis
31 Primary biliary cirrhosis.gwascatalog.bed.cut.sort.uniq
Phenotype: BMI
92 BMI.gwascatalog.bed.cut.sort.uniq
Phenotype: Type 1 diabetes
98 Type 1 diabetes.gwascatalog.bed.cut.sort.uniq
Phenotype: Longevity
47 Longevity.gwascatalog.bed.cut.sort.uniq
Phenotype: Waist-hip ratio
14 Waist-hip ratio.gwascatalog.bed.cut.sort.uniq
Phenotype: Inflammatory biomarkers
29 Inflammatory biomarkers.gwascatalog.bed.cut.sort.uniq
Phenotype: Lipid metabolism phenotypes
41 Lipid metabolism phenotypes.gwascatalog.bed.cut.sort.uniq
Phenotype: Kawasaki disease
20 Kawasaki disease.gwascatalog.bed.cut.sort.uniq
Phenotype: Pancreatic cancer
44 Pancreatic cancer.gwascatalog.bed.cut.sort.uniq
Phenotype: C-reactive protein
60 C-reactive protein.gwascatalog.bed.cut.sort.uniq
Phenotype: Celiac disease
78 Celiac disease.gwascatalog.bed.cut.sort.uniq
Phenotype: Lupus
148 Lupus.gwascatalog.bed.cut.sort.uniq
Phenotype: Myocardial infarction
13 Myocardial infarction.gwascatalog.bed.cut.sort.uniq
Phenotype: Insulin resistance
27 Insulin resistance.gwascatalog.bed.cut.sort.uniq
Phenotype: FEV1
9 FEV1.gwascatalog.bed.cut.sort.uniq
Phenotype: Fasting plasma glucose
11 Fasting plasma glucose.gwascatalog.bed.cut.sort.uniq
Phenotype: Menarche
76 Menarche.gwascatalog.bed.cut.sort.uniq
Phenotype: CardiogramplusC4D
52 CardiogramplusC4D.gwascatalog.bed.cut.sort.uniq
Converting Phenotype: Bipolar disorder
Converting Phenotype: Schizophrenia
Converting Phenotype: Crohn's disease
Converting Phenotype: Rheumatoid arthritis
Converting Phenotype: Bone mineral density
Converting Phenotype: Ulcerative colitis
Converting Phenotype: Coronary Heart
Converting Phenotype: Type 2 diabetes
Converting Phenotype: Alzheimer's disease
Converting Phenotype: Body mass index
Converting Phenotype: Multiple sclerosis
Converting Phenotype: Breast cancer
Converting Phenotype: Asthma
Converting Phenotype: Triglycerides
Converting Phenotype: Parkinson's disease
Converting Phenotype: LDL cholesterol
Converting Phenotype: Blood pressure
Converting Phenotype: Coronary Artery
Converting Phenotype: Coronary artery calcification
Converting Phenotype: Prostate cancer
Converting Phenotype: Liver enzyme levels
Converting Phenotype: Diastolic blood pressure
Converting Phenotype: Hypertension
Converting Phenotype: HDL cholesterol
Converting Phenotype: QRS duration
Converting Phenotype: Primary biliary cirrhosis
Converting Phenotype: BMI
Converting Phenotype: Type 1 diabetes
Converting Phenotype: Longevity
Converting Phenotype: Waist-hip ratio
Converting Phenotype: Inflammatory biomarkers
Converting Phenotype: Lipid metabolism phenotypes
Converting Phenotype: Kawasaki disease
Converting Phenotype: Pancreatic cancer
Converting Phenotype: C-reactive protein
Converting Phenotype: Celiac disease
Converting Phenotype: Lupus
Converting Phenotype: Myocardial infarction
Converting Phenotype: Insulin resistance
Converting Phenotype: FEV1
Converting Phenotype: Fasting plasma glucose
Converting Phenotype: Menarche
Converting Phenotype: CardiogramplusC4D
GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed
Input phenotypes:
Bipolar disorder Schizophrenia Crohn's disease Rheumatoid arthritis Bone mineral density Ulcerative colitis Coronary Heart Type 2 diabetes Alzheimer's disease Body mass index Multiple sclerosis Breast cancer Asthma Triglycerides Parkinson's disease LDL cholesterol Blood pressure Coronary Artery Coronary artery calcification Prostate cancer Liver enzyme levels Diastolic blood pressure Hypertension HDL cholesterol QRS duration Primary biliary cirrhosis BMI Type 1 diabetes Longevity Waist-hip ratio Inflammatory biomarkers Lipid metabolism phenotypes Kawasaki disease Pancreatic cancer C-reactive protein Celiac disease Lupus Myocardial infarction Insulin resistance FEV1 Fasting plasma glucose Menarche CardiogramplusC4D
Overlapping Phenotype SNPs with input bed: Bipolar disorder
Overlapping Phenotype SNPs with input bed: Schizophrenia
Overlapping Phenotype SNPs with input bed: Crohn's disease
Overlapping Phenotype SNPs with input bed: Rheumatoid arthritis
Overlapping Phenotype SNPs with input bed: Bone mineral density
Overlapping Phenotype SNPs with input bed: Ulcerative colitis
Overlapping Phenotype SNPs with input bed: Coronary Heart
Overlapping Phenotype SNPs with input bed: Type 2 diabetes
Overlapping Phenotype SNPs with input bed: Alzheimer's disease
Overlapping Phenotype SNPs with input bed: Body mass index
Overlapping Phenotype SNPs with input bed: Multiple sclerosis
Overlapping Phenotype SNPs with input bed: Breast cancer
Overlapping Phenotype SNPs with input bed: Asthma
Overlapping Phenotype SNPs with input bed: Triglycerides
Overlapping Phenotype SNPs with input bed: Parkinson's disease
Overlapping Phenotype SNPs with input bed: LDL cholesterol
Overlapping Phenotype SNPs with input bed: Blood pressure
Overlapping Phenotype SNPs with input bed: Coronary Artery
Overlapping Phenotype SNPs with input bed: Coronary artery calcification
Overlapping Phenotype SNPs with input bed: Prostate cancer
Overlapping Phenotype SNPs with input bed: Liver enzyme levels
Overlapping Phenotype SNPs with input bed: Diastolic blood pressure
Overlapping Phenotype SNPs with input bed: Hypertension
Overlapping Phenotype SNPs with input bed: HDL cholesterol
Overlapping Phenotype SNPs with input bed: QRS duration
Overlapping Phenotype SNPs with input bed: Primary biliary cirrhosis
Overlapping Phenotype SNPs with input bed: BMI
Overlapping Phenotype SNPs with input bed: Type 1 diabetes
Overlapping Phenotype SNPs with input bed: Longevity
Overlapping Phenotype SNPs with input bed: Waist-hip ratio
Overlapping Phenotype SNPs with input bed: Inflammatory biomarkers
Overlapping Phenotype SNPs with input bed: Lipid metabolism phenotypes
Overlapping Phenotype SNPs with input bed: Kawasaki disease
Overlapping Phenotype SNPs with input bed: Pancreatic cancer
Overlapping Phenotype SNPs with input bed: C-reactive protein
Overlapping Phenotype SNPs with input bed: Celiac disease
Overlapping Phenotype SNPs with input bed: Lupus
Overlapping Phenotype SNPs with input bed: Myocardial infarction
Overlapping Phenotype SNPs with input bed: Insulin resistance
Overlapping Phenotype SNPs with input bed: FEV1
Overlapping Phenotype SNPs with input bed: Fasting plasma glucose
Overlapping Phenotype SNPs with input bed: Menarche
Overlapping Phenotype SNPs with input bed: CardiogramplusC4D
Number of Overlapping Phenotype SNPs with input bed: Bipolar disorder
8 Bipolar disorder.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Schizophrenia
9 Schizophrenia.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Crohn's disease
11 Crohn's disease.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Rheumatoid arthritis
8 Rheumatoid arthritis.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Bone mineral density
6 Bone mineral density.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Ulcerative colitis
9 Ulcerative colitis.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Coronary Heart
8 Coronary Heart.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Type 2 diabetes
8 Type 2 diabetes.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Alzheimer's disease
4 Alzheimer's disease.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Body mass index
9 Body mass index.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Multiple sclerosis
10 Multiple sclerosis.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Breast cancer
7 Breast cancer.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Asthma
8 Asthma.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Triglycerides
3 Triglycerides.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Parkinson's disease
2 Parkinson's disease.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: LDL cholesterol
4 LDL cholesterol.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Blood pressure
11 Blood pressure.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Coronary Artery
1 Coronary Artery.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Coronary artery calcification
0 Coronary artery calcification.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Prostate cancer
4 Prostate cancer.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Liver enzyme levels
2 Liver enzyme levels.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Diastolic blood pressure
6 Diastolic blood pressure.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Hypertension
1 Hypertension.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: HDL cholesterol
7 HDL cholesterol.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: QRS duration
2 QRS duration.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Primary biliary cirrhosis
1 Primary biliary cirrhosis.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: BMI
1 BMI.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Type 1 diabetes
6 Type 1 diabetes.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Longevity
0 Longevity.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Waist-hip ratio
0 Waist-hip ratio.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Inflammatory biomarkers
0 Inflammatory biomarkers.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Lipid metabolism phenotypes
2 Lipid metabolism phenotypes.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Kawasaki disease
1 Kawasaki disease.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Pancreatic cancer
1 Pancreatic cancer.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: C-reactive protein
1 C-reactive protein.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Celiac disease
2 Celiac disease.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Lupus
4 Lupus.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Myocardial infarction
1 Myocardial infarction.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Insulin resistance
1 Insulin resistance.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: FEV1
0 FEV1.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Fasting plasma glucose
1 Fasting plasma glucose.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Menarche
2 Menarche.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: CardiogramplusC4D
3 CardiogramplusC4D.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Bipolar disorder
8
Number of Overlapping Phenotype SNPs with input bed: Schizophrenia
9
Number of Overlapping Phenotype SNPs with input bed: Crohn's disease
11
Number of Overlapping Phenotype SNPs with input bed: Rheumatoid arthritis
8
Number of Overlapping Phenotype SNPs with input bed: Bone mineral density
6
Number of Overlapping Phenotype SNPs with input bed: Ulcerative colitis
9
Number of Overlapping Phenotype SNPs with input bed: Coronary Heart
8
Number of Overlapping Phenotype SNPs with input bed: Type 2 diabetes
8
Number of Overlapping Phenotype SNPs with input bed: Alzheimer's disease
4
Number of Overlapping Phenotype SNPs with input bed: Body mass index
9
Number of Overlapping Phenotype SNPs with input bed: Multiple sclerosis
10
Number of Overlapping Phenotype SNPs with input bed: Breast cancer
7
Number of Overlapping Phenotype SNPs with input bed: Asthma
8
Number of Overlapping Phenotype SNPs with input bed: Triglycerides
3
Number of Overlapping Phenotype SNPs with input bed: Parkinson's disease
2
Number of Overlapping Phenotype SNPs with input bed: LDL cholesterol
4
Number of Overlapping Phenotype SNPs with input bed: Blood pressure
11
Number of Overlapping Phenotype SNPs with input bed: Coronary Artery
1
Number of Overlapping Phenotype SNPs with input bed: Coronary artery calcification
0
Number of Overlapping Phenotype SNPs with input bed: Prostate cancer
4
Number of Overlapping Phenotype SNPs with input bed: Liver enzyme levels
2
Number of Overlapping Phenotype SNPs with input bed: Diastolic blood pressure
6
Number of Overlapping Phenotype SNPs with input bed: Hypertension
1
Number of Overlapping Phenotype SNPs with input bed: HDL cholesterol
7
Number of Overlapping Phenotype SNPs with input bed: QRS duration
2
Number of Overlapping Phenotype SNPs with input bed: Primary biliary cirrhosis
1
Number of Overlapping Phenotype SNPs with input bed: BMI
1
Number of Overlapping Phenotype SNPs with input bed: Type 1 diabetes
6
Number of Overlapping Phenotype SNPs with input bed: Longevity
0
Number of Overlapping Phenotype SNPs with input bed: Waist-hip ratio
0
Number of Overlapping Phenotype SNPs with input bed: Inflammatory biomarkers
0
Number of Overlapping Phenotype SNPs with input bed: Lipid metabolism phenotypes
2
Number of Overlapping Phenotype SNPs with input bed: Kawasaki disease
1
Number of Overlapping Phenotype SNPs with input bed: Pancreatic cancer
1
Number of Overlapping Phenotype SNPs with input bed: C-reactive protein
1
Number of Overlapping Phenotype SNPs with input bed: Celiac disease
2
Number of Overlapping Phenotype SNPs with input bed: Lupus
4
Number of Overlapping Phenotype SNPs with input bed: Myocardial infarction
1
Number of Overlapping Phenotype SNPs with input bed: Insulin resistance
1
Number of Overlapping Phenotype SNPs with input bed: FEV1
0
Number of Overlapping Phenotype SNPs with input bed: Fasting plasma glucose
1
Number of Overlapping Phenotype SNPs with input bed: Menarche
2
Number of Overlapping Phenotype SNPs with input bed: CardiogramplusC4D
3
--2016-05-02 00:52:46--  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
Resolving genome.ucsc.edu (genome.ucsc.edu)... 128.114.119.132, 128.114.119.136, 128.114.119.134, ...
Connecting to genome.ucsc.edu (genome.ucsc.edu)|128.114.119.132|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1971 (1.9K) [text/plain]
Saving to: `hg19.chrom.sizes'

100%[======================================================================================================================================================>] 1,971       --.-K/s   in 0s      

2016-05-02 00:52:46 (53.9 MB/s) - `hg19.chrom.sizes' saved [1971/1971]

Human Genome size version hg19: 3137161264
99885 GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed
Coverage of BED file 85582140
Fraction of hg19 0.0272801
Overlapping Phenotype SNPs with input bed to calculate coverage: Bipolar disorder
Overlapping Phenotype SNPs with input bed to calculate coverage: Schizophrenia
Overlapping Phenotype SNPs with input bed to calculate coverage: Crohn's disease
Overlapping Phenotype SNPs with input bed to calculate coverage: Rheumatoid arthritis
Overlapping Phenotype SNPs with input bed to calculate coverage: Bone mineral density
Overlapping Phenotype SNPs with input bed to calculate coverage: Ulcerative colitis
Overlapping Phenotype SNPs with input bed to calculate coverage: Coronary Heart
Overlapping Phenotype SNPs with input bed to calculate coverage: Type 2 diabetes
Overlapping Phenotype SNPs with input bed to calculate coverage: Alzheimer's disease
Overlapping Phenotype SNPs with input bed to calculate coverage: Body mass index
Overlapping Phenotype SNPs with input bed to calculate coverage: Multiple sclerosis
Overlapping Phenotype SNPs with input bed to calculate coverage: Breast cancer
Overlapping Phenotype SNPs with input bed to calculate coverage: Asthma
Overlapping Phenotype SNPs with input bed to calculate coverage: Triglycerides
Overlapping Phenotype SNPs with input bed to calculate coverage: Parkinson's disease
Overlapping Phenotype SNPs with input bed to calculate coverage: LDL cholesterol
Overlapping Phenotype SNPs with input bed to calculate coverage: Blood pressure
Overlapping Phenotype SNPs with input bed to calculate coverage: Coronary Artery
Overlapping Phenotype SNPs with input bed to calculate coverage: Coronary artery calcification
Overlapping Phenotype SNPs with input bed to calculate coverage: Prostate cancer
Overlapping Phenotype SNPs with input bed to calculate coverage: Liver enzyme levels
Overlapping Phenotype SNPs with input bed to calculate coverage: Diastolic blood pressure
Overlapping Phenotype SNPs with input bed to calculate coverage: Hypertension
Overlapping Phenotype SNPs with input bed to calculate coverage: HDL cholesterol
Overlapping Phenotype SNPs with input bed to calculate coverage: QRS duration
Overlapping Phenotype SNPs with input bed to calculate coverage: Primary biliary cirrhosis
Overlapping Phenotype SNPs with input bed to calculate coverage: BMI
Overlapping Phenotype SNPs with input bed to calculate coverage: Type 1 diabetes
Overlapping Phenotype SNPs with input bed to calculate coverage: Longevity
Overlapping Phenotype SNPs with input bed to calculate coverage: Waist-hip ratio
Overlapping Phenotype SNPs with input bed to calculate coverage: Inflammatory biomarkers
Overlapping Phenotype SNPs with input bed to calculate coverage: Lipid metabolism phenotypes
Overlapping Phenotype SNPs with input bed to calculate coverage: Kawasaki disease
Overlapping Phenotype SNPs with input bed to calculate coverage: Pancreatic cancer
Overlapping Phenotype SNPs with input bed to calculate coverage: C-reactive protein
Overlapping Phenotype SNPs with input bed to calculate coverage: Celiac disease
Overlapping Phenotype SNPs with input bed to calculate coverage: Lupus
Overlapping Phenotype SNPs with input bed to calculate coverage: Myocardial infarction
Overlapping Phenotype SNPs with input bed to calculate coverage: Insulin resistance
Overlapping Phenotype SNPs with input bed to calculate coverage: FEV1
Overlapping Phenotype SNPs with input bed to calculate coverage: Fasting plasma glucose
Overlapping Phenotype SNPs with input bed to calculate coverage: Menarche
Overlapping Phenotype SNPs with input bed to calculate coverage: CardiogramplusC4D
Input terms no spaces: Bipolardisorder
Coverage of the input bed file that overlaps GWAS category: 12067
Fraction of hg19 3.84647e-06
Fold change: Bipolardisorder
overlap: 8
total: 444
fold:
1.8018
Input terms no spaces: Schizophrenia
Coverage of the input bed file that overlaps GWAS category: 17166
Fraction of hg19 5.47183e-06
Fold change: Schizophrenia
overlap: 9
total: 413
fold:
2.17918
Input terms no spaces: Crohnsdisease
Coverage of the input bed file that overlaps GWAS category: 16864
Fraction of hg19 5.37556e-06
Fold change: Crohnsdisease
overlap: 11
total: 228
fold:
4.82456
Input terms no spaces: Rheumatoidarthritis
Coverage of the input bed file that overlaps GWAS category: 11467
Fraction of hg19 3.65522e-06
Fold change: Rheumatoidarthritis
overlap: 8
total: 369
fold:
2.16802
Input terms no spaces: Bonemineraldensity
Coverage of the input bed file that overlaps GWAS category: 9069
Fraction of hg19 2.89083e-06
Fold change: Bonemineraldensity
overlap: 6
total: 245
fold:
2.44898
Input terms no spaces: Ulcerativecolitis
Coverage of the input bed file that overlaps GWAS category: 11016
Fraction of hg19 3.51145e-06
Fold change: Ulcerativecolitis
overlap: 9
total: 152
fold:
5.92105
Input terms no spaces: CoronaryHeart
Coverage of the input bed file that overlaps GWAS category: 10942
Fraction of hg19 3.48787e-06
Fold change: CoronaryHeart
overlap: 8
total: 143
fold:
5.59441
Input terms no spaces: Type2diabetes
Coverage of the input bed file that overlaps GWAS category: 11617
Fraction of hg19 3.70303e-06
Fold change: Type2diabetes
overlap: 8
total: 318
fold:
2.51572
Input terms no spaces: Alzheimersdisease
Coverage of the input bed file that overlaps GWAS category: 3446
Fraction of hg19 1.09845e-06
Fold change: Alzheimersdisease
overlap: 4
total: 226
fold:
1.76991
Input terms no spaces: Bodymassindex
Coverage of the input bed file that overlaps GWAS category: 13941
Fraction of hg19 4.44383e-06
Fold change: Bodymassindex
overlap: 9
total: 265
fold:
3.39623
Input terms no spaces: Multiplesclerosis
Coverage of the input bed file that overlaps GWAS category: 13565
Fraction of hg19 4.32397e-06
Fold change: Multiplesclerosis
overlap: 10
total: 215
fold:
4.65116
Input terms no spaces: Breastcancer
Coverage of the input bed file that overlaps GWAS category: 11243
Fraction of hg19 3.58381e-06
Fold change: Breastcancer
overlap: 7
total: 253
fold:
2.7668
Input terms no spaces: Asthma
Coverage of the input bed file that overlaps GWAS category: 10642
Fraction of hg19 3.39224e-06
Fold change: Asthma
overlap: 8
total: 166
fold:
4.81928
Input terms no spaces: Triglycerides
Coverage of the input bed file that overlaps GWAS category: 4047
Fraction of hg19 1.29002e-06
Fold change: Triglycerides
overlap: 3
total: 180
fold:
1.66667
Input terms no spaces: Parkinsonsdisease
Coverage of the input bed file that overlaps GWAS category: 4348
Fraction of hg19 1.38597e-06
Fold change: Parkinsonsdisease
overlap: 2
total: 105
fold:
1.90476
Input terms no spaces: LDLcholesterol
Coverage of the input bed file that overlaps GWAS category: 5546
Fraction of hg19 1.76784e-06
Fold change: LDLcholesterol
overlap: 4
total: 166
fold:
2.40964
Input terms no spaces: Bloodpressure
Coverage of the input bed file that overlaps GWAS category: 17164
Fraction of hg19 5.47119e-06
Fold change: Bloodpressure
overlap: 11
total: 266
fold:
4.13534
Input terms no spaces: CoronaryArtery
Coverage of the input bed file that overlaps GWAS category: 524
Fraction of hg19 1.6703e-07
Fold change: CoronaryArtery
overlap: 1
total: 110
fold:
0.909091
Input terms no spaces: Coronaryarterycalcification
Coverage of the input bed file that overlaps GWAS category:
Fraction of hg19 0
Fold change: Coronaryarterycalcification
overlap: 0
total: 55
fold:
0
Input terms no spaces: Prostatecancer
Coverage of the input bed file that overlaps GWAS category: 8021
Fraction of hg19 2.55677e-06
Fold change: Prostatecancer
overlap: 4
total: 225
fold:
1.77778
Input terms no spaces: Liverenzymelevels
Coverage of the input bed file that overlaps GWAS category: 1273
Fraction of hg19 4.05781e-07
Fold change: Liverenzymelevels
overlap: 2
total: 74
fold:
2.7027
Input terms no spaces: Diastolicbloodpressure
Coverage of the input bed file that overlaps GWAS category: 11469
Fraction of hg19 3.65585e-06
Fold change: Diastolicbloodpressure
overlap: 6
total: 36
fold:
16.6667
Input terms no spaces: Hypertension
Coverage of the input bed file that overlaps GWAS category: 3899
Fraction of hg19 1.24284e-06
Fold change: Hypertension
overlap: 1
total: 30
fold:
3.33333
Input terms no spaces: HDLcholesterol
Coverage of the input bed file that overlaps GWAS category: 11993
Fraction of hg19 3.82288e-06
Fold change: HDLcholesterol
overlap: 7
total: 213
fold:
3.28638
Input terms no spaces: QRSduration
Coverage of the input bed file that overlaps GWAS category: 1723
Fraction of hg19 5.49223e-07
Fold change: QRSduration
overlap: 2
total: 59
fold:
3.38983
Input terms no spaces: Primarybiliarycirrhosis
Coverage of the input bed file that overlaps GWAS category: 1199
Fraction of hg19 3.82193e-07
Fold change: Primarybiliarycirrhosis
overlap: 1
total: 38
fold:
2.63158
Input terms no spaces: BMI
Coverage of the input bed file that overlaps GWAS category: 299
Fraction of hg19 9.53091e-08
Fold change: BMI
overlap: 1
total: 96
fold:
1.04167
Input terms no spaces: Type1diabetes
Coverage of the input bed file that overlaps GWAS category: 4644
Fraction of hg19 1.48032e-06
Fold change: Type1diabetes
overlap: 6
total: 133
fold:
4.51128
Input terms no spaces: Longevity
Coverage of the input bed file that overlaps GWAS category:
Fraction of hg19 0
Fold change: Longevity
overlap: 0
total: 50
fold:
0
Input terms no spaces: Waisthipratio
Coverage of the input bed file that overlaps GWAS category:
Fraction of hg19 0
Fold change: Waisthipratio
overlap: 0
total: 14
fold:
0
Input terms no spaces: Inflammatorybiomarkers
Coverage of the input bed file that overlaps GWAS category:
Fraction of hg19 0
Fold change: Inflammatorybiomarkers
overlap: 0
total: 29
fold:
0
Input terms no spaces: Lipidmetabolismphenotypes
Coverage of the input bed file that overlaps GWAS category: 3298
Fraction of hg19 1.05127e-06
Fold change: Lipidmetabolismphenotypes
overlap: 2
total: 91
fold:
2.1978
Input terms no spaces: Kawasakidisease
Coverage of the input bed file that overlaps GWAS category: 1049
Fraction of hg19 3.34379e-07
Fold change: Kawasakidisease
overlap: 1
total: 20
fold:
5
Input terms no spaces: Pancreaticcancer
Coverage of the input bed file that overlaps GWAS category: 824
Fraction of hg19 2.62658e-07
Fold change: Pancreaticcancer
overlap: 1
total: 44
fold:
2.27273
Input terms no spaces: Creactiveprotein
Coverage of the input bed file that overlaps GWAS category: 599
Fraction of hg19 1.90937e-07
Fold change: Creactiveprotein
overlap: 1
total: 72
fold:
1.38889
Input terms no spaces: Celiacdisease
Coverage of the input bed file that overlaps GWAS category: 1723
Fraction of hg19 5.49223e-07
Fold change: Celiacdisease
overlap: 2
total: 90
fold:
2.22222
Input terms no spaces: Lupus
Coverage of the input bed file that overlaps GWAS category: 6296
Fraction of hg19 2.00691e-06
Fold change: Lupus
overlap: 4
total: 174
fold:
2.29885
Input terms no spaces: Myocardialinfarction
Coverage of the input bed file that overlaps GWAS category: 524
Fraction of hg19 1.6703e-07
Fold change: Myocardialinfarction
overlap: 1
total: 13
fold:
7.69231
Input terms no spaces: Insulinresistance
Coverage of the input bed file that overlaps GWAS category: 599
Fraction of hg19 1.90937e-07
Fold change: Insulinresistance
overlap: 1
total: 29
fold:
3.44828
Input terms no spaces: FEV1
Coverage of the input bed file that overlaps GWAS category:
Fraction of hg19 0
Fold change: FEV1
overlap: 0
total: 9
fold:
0
Input terms no spaces: Fastingplasmaglucose
Coverage of the input bed file that overlaps GWAS category: 1949
Fraction of hg19 6.21262e-07
Fold change: Fastingplasmaglucose
overlap: 1
total: 17
fold:
5.88235
Input terms no spaces: Menarche
Coverage of the input bed file that overlaps GWAS category: 1573
Fraction of hg19 5.01409e-07
Fold change: Menarche
overlap: 2
total: 82
fold:
2.43902
Input terms no spaces: CardiogramplusC4D
Coverage of the input bed file that overlaps GWAS category: 3222
Fraction of hg19 1.02704e-06
Fold change: CardiogramplusC4D
overlap: 3
total: 52
fold:
5.76923
[1] "Bipolar disorder"
[1] 1.681784e-27
               V1                   V2     V3  V4
1               1                    2      3   4
2 Bipolardisorder 1.68178406331541e-27 1.8018 444
[1] "Schizophrenia"
[1] 3.871704e-30
               V1                   V2      V3  V4
1               1                    2       3   4
2 Bipolardisorder 1.68178406331541e-27  1.8018 444
3   Schizophrenia 3.87170405753683e-30 2.17918 413
[1] "Crohn's disease"
[1] 1.835995e-40
               V1                   V2      V3  V4
1               1                    2       3   4
2 Bipolardisorder 1.68178406331541e-27  1.8018 444
3   Schizophrenia 3.87170405753683e-30 2.17918 413
4   Crohnsdisease 1.83599548155219e-40 4.82456 228
[1] "Rheumatoid arthritis"
[1] 2.513306e-28
                   V1                   V2      V3  V4
1                   1                    2       3   4
2     Bipolardisorder 1.68178406331541e-27  1.8018 444
3       Schizophrenia 3.87170405753683e-30 2.17918 413
4       Crohnsdisease 1.83599548155219e-40 4.82456 228
5 Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
[1] "Bone mineral density"
[1] 1.647057e-22
                   V1                   V2      V3  V4
1                   1                    2       3   4
2     Bipolardisorder 1.68178406331541e-27  1.8018 444
3       Schizophrenia 3.87170405753683e-30 2.17918 413
4       Crohnsdisease 1.83599548155219e-40 4.82456 228
5 Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6  Bonemineraldensity 1.64705688245707e-22 2.44898 245
[1] "Ulcerative colitis"
[1] 7.606025e-36
                   V1                   V2      V3  V4
1                   1                    2       3   4
2     Bipolardisorder 1.68178406331541e-27  1.8018 444
3       Schizophrenia 3.87170405753683e-30 2.17918 413
4       Crohnsdisease 1.83599548155219e-40 4.82456 228
5 Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6  Bonemineraldensity 1.64705688245707e-22 2.44898 245
7   Ulcerativecolitis 7.60602508855168e-36 5.92105 152
[1] "Coronary Heart"
[1] 7.778255e-32
                   V1                   V2      V3  V4
1                   1                    2       3   4
2     Bipolardisorder 1.68178406331541e-27  1.8018 444
3       Schizophrenia 3.87170405753683e-30 2.17918 413
4       Crohnsdisease 1.83599548155219e-40 4.82456 228
5 Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6  Bonemineraldensity 1.64705688245707e-22 2.44898 245
7   Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8       CoronaryHeart 7.77825526823762e-32 5.59441 143
[1] "Type 2 diabetes"
[1] 8.381307e-29
                   V1                   V2      V3  V4
1                   1                    2       3   4
2     Bipolardisorder 1.68178406331541e-27  1.8018 444
3       Schizophrenia 3.87170405753683e-30 2.17918 413
4       Crohnsdisease 1.83599548155219e-40 4.82456 228
5 Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6  Bonemineraldensity 1.64705688245707e-22 2.44898 245
7   Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8       CoronaryHeart 7.77825526823762e-32 5.59441 143
9       Type2diabetes 8.38130727686471e-29 2.51572 318
[1] "Alzheimer's disease"
[1] 1.540451e-16
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
[1] "Body mass index"
[1] 1.045041e-32
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
[1] "Multiple sclerosis"
[1] 1.073483e-37
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
[1] "Breast cancer"
[1] 9.185106e-26
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
[1] "Asthma"
[1] 2.111644e-31
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
14              Asthma 2.11164425383019e-31 4.81928 166
[1] "Triglycerides"
[1] 2.051561e-12
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
14              Asthma 2.11164425383019e-31 4.81928 166
15       Triglycerides 2.05156109382356e-12 1.66667 180
[1] "Parkinson's disease"
[1] 1.048669e-08
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
14              Asthma 2.11164425383019e-31 4.81928 166
15       Triglycerides 2.05156109382356e-12 1.66667 180
16   Parkinsonsdisease 1.04866869773618e-08 1.90476 105
[1] "LDL cholesterol"
[1] 2.978926e-16
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
14              Asthma 2.11164425383019e-31 4.81928 166
15       Triglycerides 2.05156109382356e-12 1.66667 180
16   Parkinsonsdisease 1.04866869773618e-08 1.90476 105
17      LDLcholesterol  2.9789258299434e-16 2.40964 166
[1] "Blood pressure"
[1] 1.258441e-39
                    V1                   V2      V3  V4
1                    1                    2       3   4
2      Bipolardisorder 1.68178406331541e-27  1.8018 444
3        Schizophrenia 3.87170405753683e-30 2.17918 413
4        Crohnsdisease 1.83599548155219e-40 4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28 2.16802 369
6   Bonemineraldensity 1.64705688245707e-22 2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36 5.92105 152
8        CoronaryHeart 7.77825526823762e-32 5.59441 143
9        Type2diabetes 8.38130727686471e-29 2.51572 318
10   Alzheimersdisease 1.54045091455998e-16 1.76991 226
11       Bodymassindex 1.04504050031761e-32 3.39623 265
12   Multiplesclerosis 1.07348346042233e-37 4.65116 215
13        Breastcancer 9.18510643837705e-26  2.7668 253
14              Asthma 2.11164425383019e-31 4.81928 166
15       Triglycerides 2.05156109382356e-12 1.66667 180
16   Parkinsonsdisease 1.04866869773618e-08 1.90476 105
17      LDLcholesterol  2.9789258299434e-16 2.40964 166
18       Bloodpressure  1.2584414963838e-39 4.13534 266
[1] "Coronary Artery"
[1] 1.837297e-05
                    V1                   V2       V3  V4
1                    1                    2        3   4
2      Bipolardisorder 1.68178406331541e-27   1.8018 444
3        Schizophrenia 3.87170405753683e-30  2.17918 413
4        Crohnsdisease 1.83599548155219e-40  4.82456 228
5  Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6   Bonemineraldensity 1.64705688245707e-22  2.44898 245
7    Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8        CoronaryHeart 7.77825526823762e-32  5.59441 143
9        Type2diabetes 8.38130727686471e-29  2.51572 318
10   Alzheimersdisease 1.54045091455998e-16  1.76991 226
11       Bodymassindex 1.04504050031761e-32  3.39623 265
12   Multiplesclerosis 1.07348346042233e-37  4.65116 215
13        Breastcancer 9.18510643837705e-26   2.7668 253
14              Asthma 2.11164425383019e-31  4.81928 166
15       Triglycerides 2.05156109382356e-12  1.66667 180
16   Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17      LDLcholesterol  2.9789258299434e-16  2.40964 166
18       Bloodpressure  1.2584414963838e-39  4.13534 266
19      CoronaryArtery 1.83729654937565e-05 0.909091 110
[1] "Coronary artery calcification"
[1] 1
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
[1] "Prostate cancer"
[1] 4.440158e-15
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
[1] "Liver enzyme levels"
[1] 4.447289e-10
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
[1] "Diastolic blood pressure"
[1] 4.649702e-27
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
[1] "Hypertension"
[1] 3.728386e-05
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
[1] "HDL cholesterol"
[1] 4.259507e-26
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
[1] "QRS duration"
[1] 5.161e-10
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
[1] "Primary biliary cirrhosis"
[1] 1.452313e-05
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
[1] "BMI"
[1] 9.149591e-06
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
[1] "Type 1 diabetes"
[1] 7.213722e-26
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
[1] "Longevity"
[1] 1
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
[1] "Waist-hip ratio"
[1] 1
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
[1] "Inflammatory biomarkers"
[1] 1
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
[1] "Lipid metabolism phenotypes"
[1] 4.525242e-09
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
[1] "Kawasaki disease"
[1] 6.687538e-06
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
[1] "Pancreatic cancer"
[1] 1.155682e-05
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
[1] "C-reactive protein"
[1] 1.374728e-05
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
[1] "Celiac disease"
[1] 1.208033e-09
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
[1] "Lupus"
[1] 5.982354e-16
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
[1] "Myocardial infarction"
[1] 2.171386e-06
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
[1] "Insulin resistance"
[1] 5.537143e-06
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
40           Insulinresistance 5.53714339704269e-06  3.44828  29
[1] "FEV1"
[1] 1
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
40           Insulinresistance 5.53714339704269e-06  3.44828  29
41                        FEV1                    1        0   9
[1] "Fasting plasma glucose"
[1] 1.056135e-05
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
40           Insulinresistance 5.53714339704269e-06  3.44828  29
41                        FEV1                    1        0   9
42        Fastingplasmaglucose 1.05613490176086e-05  5.88235  17
[1] "Menarche"
[1] 8.349024e-10
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
40           Insulinresistance 5.53714339704269e-06  3.44828  29
41                        FEV1                    1        0   9
42        Fastingplasmaglucose 1.05613490176086e-05  5.88235  17
43                    Menarche 8.34902391232258e-10  2.43902  82
[1] "CardiogramplusC4D"
[1] 2.394046e-14
                            V1                   V2       V3  V4
1                            1                    2        3   4
2              Bipolardisorder 1.68178406331541e-27   1.8018 444
3                Schizophrenia 3.87170405753683e-30  2.17918 413
4                Crohnsdisease 1.83599548155219e-40  4.82456 228
5          Rheumatoidarthritis 2.51330594175823e-28  2.16802 369
6           Bonemineraldensity 1.64705688245707e-22  2.44898 245
7            Ulcerativecolitis 7.60602508855168e-36  5.92105 152
8                CoronaryHeart 7.77825526823762e-32  5.59441 143
9                Type2diabetes 8.38130727686471e-29  2.51572 318
10           Alzheimersdisease 1.54045091455998e-16  1.76991 226
11               Bodymassindex 1.04504050031761e-32  3.39623 265
12           Multiplesclerosis 1.07348346042233e-37  4.65116 215
13                Breastcancer 9.18510643837705e-26   2.7668 253
14                      Asthma 2.11164425383019e-31  4.81928 166
15               Triglycerides 2.05156109382356e-12  1.66667 180
16           Parkinsonsdisease 1.04866869773618e-08  1.90476 105
17              LDLcholesterol  2.9789258299434e-16  2.40964 166
18               Bloodpressure  1.2584414963838e-39  4.13534 266
19              CoronaryArtery 1.83729654937565e-05 0.909091 110
20 Coronaryarterycalcification                    1        0  55
21              Prostatecancer 4.44015757474579e-15  1.77778 225
22           Liverenzymelevels 4.44728858620703e-10   2.7027  74
23      Diastolicbloodpressure 4.64970238571524e-27  16.6667  36
24                Hypertension 3.72838561767812e-05  3.33333  30
25              HDLcholesterol 4.25950655993691e-26  3.28638 213
26                 QRSduration 5.16099984146027e-10  3.38983  59
27     Primarybiliarycirrhosis  1.4523128624899e-05  2.63158  38
28                         BMI 9.14959075589129e-06  1.04167  96
29               Type1diabetes 7.21372191475234e-26  4.51128 133
30                   Longevity                    1        0  50
31               Waisthipratio                    1        0  14
32      Inflammatorybiomarkers                    1        0  29
33   Lipidmetabolismphenotypes 4.52524205443715e-09   2.1978  91
34             Kawasakidisease 6.68753751258792e-06        5  20
35            Pancreaticcancer 1.15568214731063e-05  2.27273  44
36            Creactiveprotein 1.37472776333785e-05  1.38889  72
37               Celiacdisease 1.20803345678883e-09  2.22222  90
38                       Lupus 5.98235366937866e-16  2.29885 174
39        Myocardialinfarction 2.17138564775674e-06  7.69231  13
40           Insulinresistance 5.53714339704269e-06  3.44828  29
41                        FEV1                    1        0   9
42        Fastingplasmaglucose 1.05613490176086e-05  5.88235  17
43                    Menarche 8.34902391232258e-10  2.43902  82
44           CardiogramplusC4D  2.3940460110921e-14  5.76923  52
         V1          V2          V3          V4 
"character"   "numeric"   "numeric"   "numeric" 
                            V1           V2        V3  V4
2              Bipolardisorder 1.681784e-27  1.801800 444
3                Schizophrenia 3.871704e-30  2.179180 413
4                Crohnsdisease 1.835995e-40  4.824560 228
5          Rheumatoidarthritis 2.513306e-28  2.168020 369
6           Bonemineraldensity 1.647057e-22  2.448980 245
7            Ulcerativecolitis 7.606025e-36  5.921050 152
8                CoronaryHeart 7.778255e-32  5.594410 143
9                Type2diabetes 8.381307e-29  2.515720 318
10           Alzheimersdisease 1.540451e-16  1.769910 226
11               Bodymassindex 1.045041e-32  3.396230 265
12           Multiplesclerosis 1.073483e-37  4.651160 215
13                Breastcancer 9.185106e-26  2.766800 253
14                      Asthma 2.111644e-31  4.819280 166
15               Triglycerides 2.051561e-12  1.666670 180
16           Parkinsonsdisease 1.048669e-08  1.904760 105
17              LDLcholesterol 2.978926e-16  2.409640 166
18               Bloodpressure 1.258441e-39  4.135340 266
19              CoronaryArtery 1.837297e-05  0.909091 110
20 Coronaryarterycalcification 1.000000e+00  0.000000  55
21              Prostatecancer 4.440158e-15  1.777780 225
22           Liverenzymelevels 4.447289e-10  2.702700  74
23      Diastolicbloodpressure 4.649702e-27 16.666700  36
24                Hypertension 3.728386e-05  3.333330  30
25              HDLcholesterol 4.259507e-26  3.286380 213
26                 QRSduration 5.161000e-10  3.389830  59
27     Primarybiliarycirrhosis 1.452313e-05  2.631580  38
28                         BMI 9.149591e-06  1.041670  96
29               Type1diabetes 7.213722e-26  4.511280 133
30                   Longevity 1.000000e+00  0.000000  50
31               Waisthipratio 1.000000e+00  0.000000  14
32      Inflammatorybiomarkers 1.000000e+00  0.000000  29
33   Lipidmetabolismphenotypes 4.525242e-09  2.197800  91
34             Kawasakidisease 6.687538e-06  5.000000  20
35            Pancreaticcancer 1.155682e-05  2.272730  44
36            Creactiveprotein 1.374728e-05  1.388890  72
37               Celiacdisease 1.208033e-09  2.222220  90
38                       Lupus 5.982354e-16  2.298850 174
39        Myocardialinfarction 2.171386e-06  7.692310  13
40           Insulinresistance 5.537143e-06  3.448280  29
41                        FEV1 1.000000e+00  0.000000   9
42        Fastingplasmaglucose 1.056135e-05  5.882350  17
43                    Menarche 8.349024e-10  2.439020  82
44           CardiogramplusC4D 2.394046e-14  5.769230  52
                            V1       V2        V3  V4
2              Bipolardisorder 61.64994  1.801800 444
3                Schizophrenia 67.72386  2.179180 413
4                Crohnsdisease 91.49582  4.824560 228
5          Rheumatoidarthritis 63.55078  2.168020 369
6           Bonemineraldensity 50.15788  2.448980 245
7            Ulcerativecolitis 80.86412  5.921050 152
8                CoronaryHeart 71.63139  5.594410 143
9                Type2diabetes 64.64896  2.515720 318
10           Alzheimersdisease 36.40929  1.769910 226
11               Bodymassindex 73.63867  3.396230 265
12           Multiplesclerosis 85.12474  4.651160 215
13                Breastcancer 57.64963  2.766800 253
14                      Asthma 70.63267  4.819280 166
15               Triglycerides 26.91242  1.666670 180
16           Parkinsonsdisease 18.37316  1.904760 105
17              LDLcholesterol 35.74980  2.409640 166
18               Bloodpressure 89.57094  4.135340 266
19              CoronaryArtery 10.90463  0.909091 110
20 Coronaryarterycalcification  0.00000  0.000000  55
21              Prostatecancer 33.04809  1.777780 225
22           Liverenzymelevels 21.53356  2.702700  74
23      Diastolicbloodpressure 60.63299 16.666700  36
24                Hypertension 10.19695  3.333330  30
25              HDLcholesterol 58.41806  3.286380 213
26                 QRSduration 21.38472  3.389830  59
27     Primarybiliarycirrhosis 11.13977  2.631580  38
28                         BMI 11.60180  1.041670  96
29               Type1diabetes 57.89123  4.511280 133
30                   Longevity  0.00000  0.000000  50
31               Waisthipratio  0.00000  0.000000  14
32      Inflammatorybiomarkers  0.00000  0.000000  29
33   Lipidmetabolismphenotypes 19.21359  2.197800  91
34             Kawasakidisease 11.91526  5.000000  20
35            Pancreaticcancer 11.36823  2.272730  44
36            Creactiveprotein 11.19467  1.388890  72
37               Celiacdisease 20.53427  2.222220  90
38                       Lupus 35.05255  2.298850 174
39        Myocardialinfarction 13.04015  7.692310  13
40           Insulinresistance 12.10403  3.448280  29
41                        FEV1  0.00000  0.000000   9
42        Fastingplasmaglucose 11.45831  5.882350  17
43                    Menarche 20.90371  2.439020  82
44           CardiogramplusC4D 31.36321  5.769230  52
                            V1       V2        V3  V4
2              Bipolardisorder 61.64994  1.801800 444
3                Schizophrenia 67.72386  2.179180 413
4                Crohnsdisease 91.49582  4.824560 228
5          Rheumatoidarthritis 63.55078  2.168020 369
6           Bonemineraldensity 50.15788  2.448980 245
7            Ulcerativecolitis 80.86412  5.921050 152
8                CoronaryHeart 71.63139  5.594410 143
9                Type2diabetes 64.64896  2.515720 318
10           Alzheimersdisease 36.40929  1.769910 226
11               Bodymassindex 73.63867  3.396230 265
12           Multiplesclerosis 85.12474  4.651160 215
13                Breastcancer 57.64963  2.766800 253
14                      Asthma 70.63267  4.819280 166
15               Triglycerides 26.91242  1.666670 180
16           Parkinsonsdisease 18.37316  1.904760 105
17              LDLcholesterol 35.74980  2.409640 166
18               Bloodpressure 89.57094  4.135340 266
19              CoronaryArtery 10.90463  0.909091 110
20 Coronaryarterycalcification  0.00000  0.000000  55
21              Prostatecancer 33.04809  1.777780 225
22           Liverenzymelevels 21.53356  2.702700  74
23      Diastolicbloodpressure 60.63299 16.666700  36
24                Hypertension 10.19695  3.333330  30
25              HDLcholesterol 58.41806  3.286380 213
26                 QRSduration 21.38472  3.389830  59
27     Primarybiliarycirrhosis 11.13977  2.631580  38
28                         BMI 11.60180  1.041670  96
29               Type1diabetes 57.89123  4.511280 133
30                   Longevity  0.00000  0.000000  50
31               Waisthipratio  0.00000  0.000000  14
32      Inflammatorybiomarkers  0.00000  0.000000  29
33   Lipidmetabolismphenotypes 19.21359  2.197800  91
34             Kawasakidisease 11.91526  5.000000  20
35            Pancreaticcancer 11.36823  2.272730  44
36            Creactiveprotein 11.19467  1.388890  72
37               Celiacdisease 20.53427  2.222220  90
38                       Lupus 35.05255  2.298850 174
39        Myocardialinfarction 13.04015  7.692310  13
40           Insulinresistance 12.10403  3.448280  29
41                        FEV1  0.00000  0.000000   9
42        Fastingplasmaglucose 11.45831  5.882350  17
43                    Menarche 20.90371  2.439020  82
44           CardiogramplusC4D 31.36321  5.769230  52
                                  V2        V3  V4                   V5
Bipolardisorder             61.64994  1.801800 444                Brain
Schizophrenia               67.72386  2.179180 413                Brain
Crohnsdisease               91.49582  4.824560 228 Chronic Inflammatory
Rheumatoidarthritis         63.55078  2.168020 369 Chronic Inflammatory
Bonemineraldensity          50.15788  2.448980 245                Other
Ulcerativecolitis           80.86412  5.921050 152 Chronic Inflammatory
CoronaryHeart               71.63139  5.594410 143       Cardiovascular
Type2diabetes               64.64896  2.515720 318                Other
Alzheimersdisease           36.40929  1.769910 226                Brain
Bodymassindex               73.63867  3.396230 265                Other
Multiplesclerosis           85.12474  4.651160 215                Brain
Breastcancer                57.64963  2.766800 253               Cancer
Asthma                      70.63267  4.819280 166                Other
Triglycerides               26.91242  1.666670 180                Other
Parkinsonsdisease           18.37316  1.904760 105                Brain
LDLcholesterol              35.74980  2.409640 166                Other
Bloodpressure               89.57094  4.135340 266                Other
CoronaryArtery              10.90463  0.909091 110       Cardiovascular
Coronaryarterycalcification  0.00000  0.000000  55       Cardiovascular
Prostatecancer              33.04809  1.777780 225               Cancer
Liverenzymelevels           21.53356  2.702700  74                Other
Diastolicbloodpressure      60.63299 16.666700  36                Other
Hypertension                10.19695  3.333330  30                Other
HDLcholesterol              58.41806  3.286380 213                Other
QRSduration                 21.38472  3.389830  59                Other
Primarybiliarycirrhosis     11.13977  2.631580  38                Other
BMI                         11.60180  1.041670  96                Other
Type1diabetes               57.89123  4.511280 133                Other
Longevity                    0.00000  0.000000  50                Other
Waisthipratio                0.00000  0.000000  14                Other
Inflammatorybiomarkers       0.00000  0.000000  29                Other
Lipidmetabolismphenotypes   19.21359  2.197800  91                Other
Kawasakidisease             11.91526  5.000000  20                Other
Pancreaticcancer            11.36823  2.272730  44               Cancer
Creactiveprotein            11.19467  1.388890  72                Other
Celiacdisease               20.53427  2.222220  90                Other
Lupus                       35.05255  2.298850 174 Chronic Inflammatory
Myocardialinfarction        13.04015  7.692310  13       Cardiovascular
Insulinresistance           12.10403  3.448280  29                Other
FEV1                         0.00000  0.000000   9                Other
Fastingplasmaglucose        11.45831  5.882350  17                Other
Menarche                    20.90371  2.439020  82                Other
CardiogramplusC4D           31.36321  5.769230  52       Cardiovascular
            Bipolardisorder               Schizophrenia 
                       TRUE                        TRUE 
              Crohnsdisease         Rheumatoidarthritis 
                       TRUE                        TRUE 
         Bonemineraldensity           Ulcerativecolitis 
                       TRUE                        TRUE 
              CoronaryHeart               Type2diabetes 
                       TRUE                        TRUE 
          Alzheimersdisease               Bodymassindex 
                       TRUE                        TRUE 
          Multiplesclerosis                Breastcancer 
                       TRUE                        TRUE 
                     Asthma               Triglycerides 
                       TRUE                        TRUE 
          Parkinsonsdisease              LDLcholesterol 
                       TRUE                        TRUE 
              Bloodpressure              CoronaryArtery 
                       TRUE                        TRUE 
Coronaryarterycalcification              Prostatecancer 
                      FALSE                        TRUE 
          Liverenzymelevels      Diastolicbloodpressure 
                       TRUE                        TRUE 
               Hypertension              HDLcholesterol 
                       TRUE                        TRUE 
                QRSduration     Primarybiliarycirrhosis 
                       TRUE                        TRUE 
                        BMI               Type1diabetes 
                       TRUE                        TRUE 
                  Longevity               Waisthipratio 
                      FALSE                       FALSE 
     Inflammatorybiomarkers   Lipidmetabolismphenotypes 
                      FALSE                        TRUE 
            Kawasakidisease            Pancreaticcancer 
                       TRUE                        TRUE 
           Creactiveprotein               Celiacdisease 
                       TRUE                        TRUE 
                      Lupus        Myocardialinfarction 
                       TRUE                        TRUE 
          Insulinresistance                        FEV1 
                       TRUE                       FALSE 
       Fastingplasmaglucose                    Menarche 
                       TRUE                        TRUE 
          CardiogramplusC4D 
                       TRUE 
                                V2        V3  V4                   V5
Bipolardisorder           61.64994  1.801800 444                Brain
Schizophrenia             67.72386  2.179180 413                Brain
Crohnsdisease             91.49582  4.824560 228 Chronic Inflammatory
Rheumatoidarthritis       63.55078  2.168020 369 Chronic Inflammatory
Bonemineraldensity        50.15788  2.448980 245                Other
Ulcerativecolitis         80.86412  5.921050 152 Chronic Inflammatory
CoronaryHeart             71.63139  5.594410 143       Cardiovascular
Type2diabetes             64.64896  2.515720 318                Other
Alzheimersdisease         36.40929  1.769910 226                Brain
Bodymassindex             73.63867  3.396230 265                Other
Multiplesclerosis         85.12474  4.651160 215                Brain
Breastcancer              57.64963  2.766800 253               Cancer
Asthma                    70.63267  4.819280 166                Other
Triglycerides             26.91242  1.666670 180                Other
Parkinsonsdisease         18.37316  1.904760 105                Brain
LDLcholesterol            35.74980  2.409640 166                Other
Bloodpressure             89.57094  4.135340 266                Other
CoronaryArtery            10.90463  0.909091 110       Cardiovascular
Prostatecancer            33.04809  1.777780 225               Cancer
Liverenzymelevels         21.53356  2.702700  74                Other
Diastolicbloodpressure    60.63299 16.666700  36                Other
Hypertension              10.19695  3.333330  30                Other
HDLcholesterol            58.41806  3.286380 213                Other
QRSduration               21.38472  3.389830  59                Other
Primarybiliarycirrhosis   11.13977  2.631580  38                Other
BMI                       11.60180  1.041670  96                Other
Type1diabetes             57.89123  4.511280 133                Other
Lipidmetabolismphenotypes 19.21359  2.197800  91                Other
Kawasakidisease           11.91526  5.000000  20                Other
Pancreaticcancer          11.36823  2.272730  44               Cancer
Creactiveprotein          11.19467  1.388890  72                Other
Celiacdisease             20.53427  2.222220  90                Other
Lupus                     35.05255  2.298850 174 Chronic Inflammatory
Myocardialinfarction      13.04015  7.692310  13       Cardiovascular
Insulinresistance         12.10403  3.448280  29                Other
Fastingplasmaglucose      11.45831  5.882350  17                Other
Menarche                  20.90371  2.439020  82                Other
CardiogramplusC4D         31.36321  5.769230  52       Cardiovascular
                              LogP        FC Phenotype SNPs
Bipolardisorder           61.64994  1.801800            444
Schizophrenia             67.72386  2.179180            413
Crohnsdisease             91.49582  4.824560            228
Rheumatoidarthritis       63.55078  2.168020            369
Bonemineraldensity        50.15788  2.448980            245
Ulcerativecolitis         80.86412  5.921050            152
CoronaryHeart             71.63139  5.594410            143
Type2diabetes             64.64896  2.515720            318
Alzheimersdisease         36.40929  1.769910            226
Bodymassindex             73.63867  3.396230            265
Multiplesclerosis         85.12474  4.651160            215
Breastcancer              57.64963  2.766800            253
Asthma                    70.63267  4.819280            166
Triglycerides             26.91242  1.666670            180
Parkinsonsdisease         18.37316  1.904760            105
LDLcholesterol            35.74980  2.409640            166
Bloodpressure             89.57094  4.135340            266
CoronaryArtery            10.90463  0.909091            110
Prostatecancer            33.04809  1.777780            225
Liverenzymelevels         21.53356  2.702700             74
Diastolicbloodpressure    60.63299 16.666700             36
Hypertension              10.19695  3.333330             30
HDLcholesterol            58.41806  3.286380            213
QRSduration               21.38472  3.389830             59
Primarybiliarycirrhosis   11.13977  2.631580             38
BMI                       11.60180  1.041670             96
Type1diabetes             57.89123  4.511280            133
Lipidmetabolismphenotypes 19.21359  2.197800             91
Kawasakidisease           11.91526  5.000000             20
Pancreaticcancer          11.36823  2.272730             44
Creactiveprotein          11.19467  1.388890             72
Celiacdisease             20.53427  2.222220             90
Lupus                     35.05255  2.298850            174
Myocardialinfarction      13.04015  7.692310             13
Insulinresistance         12.10403  3.448280             29
Fastingplasmaglucose      11.45831  5.882350             17
Menarche                  20.90371  2.439020             82
CardiogramplusC4D         31.36321  5.769230             52
                                      Category
Bipolardisorder                          Brain
Schizophrenia                            Brain
Crohnsdisease             Chronic Inflammatory
Rheumatoidarthritis       Chronic Inflammatory
Bonemineraldensity                       Other
Ulcerativecolitis         Chronic Inflammatory
CoronaryHeart                   Cardiovascular
Type2diabetes                            Other
Alzheimersdisease                        Brain
Bodymassindex                            Other
Multiplesclerosis                        Brain
Breastcancer                            Cancer
Asthma                                   Other
Triglycerides                            Other
Parkinsonsdisease                        Brain
LDLcholesterol                           Other
Bloodpressure                            Other
CoronaryArtery                  Cardiovascular
Prostatecancer                          Cancer
Liverenzymelevels                        Other
Diastolicbloodpressure                   Other
Hypertension                             Other
HDLcholesterol                           Other
QRSduration                              Other
Primarybiliarycirrhosis                  Other
BMI                                      Other
Type1diabetes                            Other
Lipidmetabolismphenotypes                Other
Kawasakidisease                          Other
Pancreaticcancer                        Cancer
Creactiveprotein                         Other
Celiacdisease                            Other
Lupus                     Chronic Inflammatory
Myocardialinfarction            Cardiovascular
Insulinresistance                        Other
Fastingplasmaglucose                     Other
Menarche                                 Other
CardiogramplusC4D               Cardiovascular
Loading required package: methods
Loading required package: grid
Loading required package: quadprog
Loading required package: proto
Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
null device 
          1 
