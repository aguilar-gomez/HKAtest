#!/usr/bin/env python3
# Author : Diana Aguilar
'''
HKA test for 3 populations: Calculate pvalues of chi-squared test comparing the polymorphic and fixed sites of a user defined window vs genome-wide
Output: chromosome	window_start	window_end	chi2_pop1	chi2_pop2	chi2_pop3	
'''
import pandas as pd
import scipy.stats as ss
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fixed_af", type=float, help="frequency of allele to be considered fixed, for low coverage data it is not recommended to use 1")
parser.add_argument("pop1", type=str, help="name of population 1 (prefix of maf file)")
parser.add_argument("pop2", type=str, help="name of population 2 (prefix of maf file)")
parser.add_argument("pop3", type=str, help="name of population 3 (prefix of maf file)")
parser.add_argument("window_size", type=int, help="Size of base pair window")
parser.add_argument("slide", type=int, help="Size of step for sliding")
args = parser.parse_args()

def fixed(af1,af2,af3):
	'''
	Returns 1 if the allele is fixed in pop1, 0 otherwise
	'''
	af1_1 = (af1>fixed_af) & (af2<1-fixed_af) & (af3<1-fixed_af)
	af1_0 = (af1<1-fixed_af) & (af2>fixed_af) & (af3>fixed_af)
	af1_fixed = af1_1 | af1_0
	return af1_fixed.astype(int)

def poly(af):
	'''
	Returns 1 if the allele is polymorphic in pop1, 0 otherwise
	'''
	af_p = (af>=1-fixed_af) & (af<=fixed_af)
	return af_p.astype(int)

def genome_wide(pop):
	'''
	Calculates the genome-wide proportion of polymorphic and fixed sites
	'''
	pol = sum( p1p2p3[pop+"_seg"])
	fix = sum( p1p2p3[pop+"_fix"])
	pol_proportion = pol /(pol + fix)
	fix_proportion = fix /(pol + fix)
	return 	pol_proportion,fix_proportion

def chi_test_w(counts_w,pop,pol_prop,fix_prop):
	'''
	Performs a chi-squared test (goodness of fit) of a window vs the expected value based on the genome wide distribution
	'''
	n_seg = counts_w[pop+"_seg"]
	n_fix = counts_w[pop+"_fix"]
	n_sites = n_seg + n_fix
	return ss.chisquare([n_seg, n_fix],[n_sites*pol_prop,n_sites*fix_prop])

outfile = args.pop1 + "." + args.pop2 + "." + args.pop3 + "." + "HKA_test"
fixed_af = args.fixed_af
pop1=pd.read_csv(args.pop1 + ".mafs",sep="\t")
pop2=pd.read_csv(args.pop2 + ".mafs",sep="\t")
pop3=pd.read_csv(args.pop3 + ".mafs",sep="\t")

pop1pop2=pop1.merge(pop2,on=["chromo","position","major","minor","ref"],suffixes=["pop1","pop2"])
p1p2p3=pop1pop2.merge(pop3,on=["chromo","position","major","minor","ref"])

print("Analyzing "+str(len(p1p2p3))+" sites found in the 3 populations")

#Fixed D
p1p2p3["pop1_fix"]=fixed(p1p2p3["knownEMpop1"],p1p2p3["knownEMpop2"],p1p2p3["knownEM"])
p1p2p3["pop2_fix"]=fixed(p1p2p3["knownEMpop2"],p1p2p3["knownEMpop1"],p1p2p3["knownEM"])
p1p2p3["pop3_fix"]=fixed(p1p2p3["knownEM"],p1p2p3["knownEMpop1"],p1p2p3["knownEMpop2"])

#Segregating/polymorphic P
p1p2p3["pop1_seg"]=poly(p1p2p3["knownEMpop1"])
p1p2p3["pop2_seg"]=poly(p1p2p3["knownEMpop2"])
p1p2p3["pop3_seg"]=poly(p1p2p3["knownEM"])

#Genome wide proportions
p1_pol_prop, p1_fix_prop = genome_wide("pop1")
p2_pol_prop, p2_fix_prop = genome_wide("pop2")
p3_pol_prop, p3_fix_prop = genome_wide("pop3")

print("Proportion of polymorphic and fixed "+args.pop1+":",p1_pol_prop, p1_fix_prop)
print("Proportion of polymorphic and fixed "+args.pop2+":",p2_pol_prop, p2_fix_prop)
print("Proportion of polymorphic and fixed "+args.pop3+":",p3_pol_prop, p3_fix_prop)

#Make windows and do HKA test each window
chr_sizes = p1p2p3.groupby(["chromo"]).count()

features=["pop1_fix","pop1_seg","pop2_fix","pop2_seg","pop3_fix","pop3_seg"]
with open(outfile,'w') as f:
	f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format("chromosome", "w_start", "w_end","chi2_p1","chi2_2","chi_p3"))
	for chromosome in chr_sizes.index:
		smalldf = p1p2p3[p1p2p3["chromo"]==chromosome].reset_index()
		w_start = 1 #angsd is one-based
		while w_start+args.window_size <= smalldf["position"].max():
			w_end = w_start+args.window_size
			counts_w = smalldf[(smalldf["position"] >= w_start) & (smalldf["position"] <= w_end)][features].sum()
			chi_p1 = chi_test_w(counts_w,"pop1",p1_pol_prop, p1_fix_prop)
			chi_p2 = chi_test_w(counts_w,"pop2",p2_pol_prop, p2_fix_prop)
			chi_p3 = chi_test_w(counts_w,"pop3",p3_pol_prop, p3_fix_prop)
			f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(chromosome, w_start, w_end, chi_p1[1], chi_p2[1], chi_p3[1]) )
			w_start += args.slide
f.close()
