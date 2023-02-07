#!/usr/bin/env python3
# coding: utf-8
# Author : Diana Aguilar

import pandas as pd
import scipy.stats as ss
import argparse
import numpy as np

def fixed(af1,af2, minor1, minor2):
	'''
	Returns 1 if the allele is fixed in pop1, 0 otherwise
	'''
	af1_1 = (af1>fixed_af) & (af2<1-fixed_af)
	af1_0 = (af1<1-fixed_af) & (af2>fixed_af)
	af1_fixed = af1_1 | af1_0
	
	#Allele is not segregating in pop2 and the minor allele differs 
	af2_fixed = (af2<1-fixed_af) | (af2>fixed_af)
	triallelic= (!af1_fixed & af2_fixed) & (minor1 != minor2)
		
	fixed = triallelic | af1_fixed
	
	return fixed.astype(int)

def poly(af1, af2):
	'''
	Returns 1 if the allele is polymorphic in pop1, 0 otherwise
	'''
	af1_p = (af1>=1-fixed_af) & (af1<=fixed_af) 
	af2_fix = (af2<1-fixed_af) | (af2>fixed_af)
	af_p = af1_p & af2_fix
	return af_p.astype(int)

def genome_wide(pop):
	'''
	Calculates the genome-wide proportion of polymorphic and fixed sites
	'''
	pol = sum( p1p2[pop+"_seg"])
	fix = sum( p1p2[pop+"_fix"])
	pol_proportion = pol /(pol + fix)
	fix_proportion = fix /(pol + fix)
	return 	pol_proportion,fix_proportion

def genome_count(pop):
	'''
	Calculates the genome-wide count of polymorphic and fixed sites
	'''
	pol = sum( p1p2[pop+"_seg"])
	fix = sum( p1p2[pop+"_fix"])
	return 	pol,fix

def chi_contingency(counts_w,pop,pol,fix):
    '''
    Performs a chi-squared contingency test 
    '''
    n_seg = counts_w[pop+"_seg"]
    n_fix = counts_w[pop+"_fix"]
    table=np.array([[n_seg, n_fix],[pol,fix]])
    if n_seg+n_fix>0:
	#Apply Yates correction to chi2
        chi2, p, dof, ex = ss.chi2_contingency(table,correction=True)
        return [chi2,p,n_seg, n_fix,pol,fix]
    else:
        return ["error",1,1,"na"]+[n_seg, n_fix,pol,fix]

parser = argparse.ArgumentParser()
parser.add_argument("fixed_af", type=float, help="frequency of allele to be considered fixed, for low coverage data it is not recommended to use 1")
parser.add_argument("pop1", type=str, help="name of population 1 (prefix of maf file)")
parser.add_argument("pop2", type=str, help="name of population 2 (prefix of maf file)")
args = parser.parse_args()

outfile = args.pop1 + "." + args.pop2 + ".fa" + str(args.fixed_af) + ".HKAgene_test"
fixed_af = args.fixed_af
pop1=pd.read_csv(args.pop1 + ".mafs",sep="\t")
pop2=pd.read_csv(args.pop2 + ".mafs",sep="\t")

#Added minor in the fileds they have to match to secure biallelic positions
#p1p2=pop1.merge(pop2,on=["chromo","position","major","ref","minor"],suffixes=["pop1","pop2"])

p1p2=pop1.merge(pop2,on=["chromo","position","major","ref"],suffixes=["pop1","pop2"])

print("Analyzing "+str(len(p1p2))+" sites found in both populations")


#Fixed D
p1p2["pop1_fix"]=fixed(p1p2["knownEMpop1"],p1p2["knownEMpop2"],p1p2["minorpop1"],p1p2["minorpop2"])
p1p2["pop2_fix"]=fixed(p1p2["knownEMpop2"],p1p2["knownEMpop1"],p1p2["minorpop2"],p1p2["minorpop1"])

#Segregating/polymorphic P
p1p2["pop1_seg"]=poly(p1p2["knownEMpop1"],p1p2["knownEMpop2"])
p1p2["pop2_seg"]=poly(p1p2["knownEMpop2"],p1p2["knownEMpop1"])

#Genome wide proportions
p1_pol_prop, p1_fix_prop = genome_wide("pop1")
p2_pol_prop, p2_fix_prop = genome_wide("pop2")


print("Proportion of polymorphic and fixed "+args.pop1+":",p1_pol_prop, p1_fix_prop)
print("Proportion of polymorphic and fixed "+args.pop2+":",p2_pol_prop, p2_fix_prop)

#Genome wide counts
p1_pol, p1_fix = genome_count("pop1")
p2_pol, p2_fix = genome_count("pop2")
print("Count of polymorphic and fixed "+args.pop1+":",p1_pol, p1_fix)
print("Count of polymorphic and fixed "+args.pop2+":",p2_pol, p2_fix)


print("Total sites "+args.pop1+":", p1_pol+p1_fix)
print("Total sites "+args.pop2+":",p2_pol+p2_fix)

chr_sizes = p1p2.groupby(["chromo"]).count()
features=["pop1_fix","pop1_seg","pop2_fix","pop2_seg"]


f = open(outfile,"w")
columns=["gene","chi2","pvalue","poly_gene","fix_gene","poly_transcriptome","fix_transcriptome"]
f.write("\t".join(columns))



for chromosome in chr_sizes.index:
        smalldf = p1p2[p1p2["chromo"]==chromosome].reset_index()
        counts_w = smalldf[features].sum()
        if p1_pol!=0 or p1_fix!=0:
                chi_p1 = chi_contingency(counts_w,"pop1",p1_pol, p1_fix)
                if chi_p1[0]!="error":
                        f.write("\n")
                        f.write("\t".join([chromosome] + [str(x) for x in chi_p1] ) )

f.close()
