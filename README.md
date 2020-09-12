# HKAtest
The following scripts are available:
1) *HKAtest.py* works on genome assemblies and calculates HKA on windows of user defined size, it works on 3 populations
2) *HKAgenes.py* works on transcriptomes and calculated HKA comparing each gene to the transcriptome, each gene gets a pvalue, for 2 populations
3) *chi-correction.R* is a script that calculates the inflation factor from the chi-squared statistic, corrects the test statistics and recalculates the p-value

## Installation
git clone https://github.com/aguilar-gomez/HKAtest.git

## HKAtest.py
### Usage 
This script performs an HKA test over 3 populations. It uses the allele frequencies generated in angsd. It uses as input the output files of -doMaf from angsd (unzipped pop.maf.gz). 
Note: make sure to use the same reference for all the bams and when running angsd, use the option -doMajorMinor 4 which uses the reference as the major allele (check angsd github and wiki in case this gets updated).

Usage: HKAtest.py [-h] fixed_af pop1 pop2 pop3 window_size slide

positional arguments:
  - fixed_af     frequency of allele to be considered fixed, for low coverage data it is not recommended to use 1 (i.e. .95)
  - pop1         name of population 1 (prefix of maf file)
  - pop2         name of population 2 (prefix of maf file)
  - pop3         name of population 3 (prefix of maf file)
  - window_size  Size of base pair window
  - slide        Size of step for sliding

optional arguments:
  -h, --help   show this help message and exit


### Output 
It prints the number of analyzed sites (present in the three populations)

It also prints the observed proportions per population of fixed vs polymorphic sites

Example:

Analyzing 28779842 sites found in the 3 populations

- Proportion of polymorphic and fixed pop1: *polymorphic_proportion_pop1 fixed_proportion_pop1= 1 - polymorphic_proportion_pop1*
- Proportion of polymorphic and fixed pop2: *polymorphic_proportion_pop2 fixed_proportion_pop2*
- Proportion of polymorphic and fixed pop3: *polymorphic_proportion_pop3 fixed_proportion_pop3*

It generates a file (pop1pop2pop3.HKA_test) with the following columns:
  - chromosome
  - w_start	: window start
  - w_end : window end
  - chi2_p1: p-value of chi-squared test using expected and observed fixed differences in pop1
  - chi2_p2:  p-value of chi-squared test using expected and observed fixed differences in pop2
  - chi_p3: p-value of chi-squared test using expected and observed fixed differences in pop3
  
## HKAgenes.py
### Usage 
This script performs an HKA test with 2 populations. It uses the allele frequencies generated in angsd. It uses as input the output files of -doMaf from angsd (unzipped pop.maf.gz). 
Note: make sure to use the same reference for all the bams and when running angsd, use the option -doMajorMinor 4 which uses the reference as the major allele (check angsd github and wiki in case this gets updated).

usage: HKAgenes.py [-h] fixed_af pop1 pop2

positional arguments:
  - fixed_af    frequency of allele to be considered fixed, for low coverage data it is not recommended to use 1
  - pop1        name of population 1 (prefix of maf file)
  - pop2        name of population 2 (prefix of maf file)

optional arguments:
  -h, --help  show this help message and exit

### Output 
It prints the number of analyzed sites (present in the both populations)

It also prints the observed proportions per population of fixed vs polymorphic sites

It generates a file (pop1.pop2.fa.HKA_genes) with the following columns:
  - gene
  - chi: chi-squared contingency table test of gene vs transcritome-wide polymorphic and fixed differences in pop1
  - pvalue:  p-value of chi-squared test
  - degrees of freedom
  - number of polymorphic positions in that gene
  - number of fixed positions in that gene
  - number of polymorphic positions in all transcriptome
  - number of fixed positions in all transcriptome
  
## chi-correction.R
### Usage 
This script takes the output file of *HKAgenes.py* and outputs a corrected file with the same name +"_corrected.tab". 

This script requires the library snpStats.

usage: chi-correction.R outfile_HKAgenes name min_sites

arguments:
  - outfile_HKAgenes
  - name          name suffix for plots generated  
  - min_sites     minimum number of sites to include the gene in the analysis*
  
*Chi-squared contingency test may not be considered valid if observed and expected values do not have a minimum number of observations. Sites = observed polymorphic + fixed. Eliminating these genes will most likely delete tests that are not significant anyway, then the correction is going to be done in a distribution that has less p-values~1, the correction might be more strict. If you want to include ALL genes regardless, set this parameter to 0


### Output 
It generates a file (pop1.pop2.fa.HKA_genes_corrected.tab) with the same columns as the output from *HKAgenes.py* plus:
  - sites: observed polymorphic + observed fixed
  - corrected: corrected chi-squared test statistic
  - pval2: corrected p-value, calculated from the corrected chi-squared statistic
  
It also outputs two plots:
 - qqplot_chiHKA.name.png: expected vs observed chi-squared test statistic before correction
 - qqplot_chiHKA_corrected.name.png: expected vs observed chi-squared test statistic after correction
