# OOP implementation of the informative read depth call step for SIBSAM
# Designed to run 1 chrm arm at a time!
# Will produce the _windows files for SIBSAM but they will need to be merged (all the chrm arms into a single windows file)
# Example command: python3 ancestry_difference_sibsam.py -c X -p1 ../../parentals/EF96N_chrXtest.vcf -p2 ../../parentals/ZI192N_chrXtest.vcf -i1 ../../parentals/EF96Nrealign_INDELS_chrX.vcf -i2 ../../parentals/ZI192Nrealign_INDELS_chrX.vcf -f1 ../../C2_remappedrealign_chr4_paired_linear_test.sam -f2 ../../W2_remappedrealign_chr4_paired_linear_test.sam -w /home/tiago/Documents/PoolLab/coldSIBSAM/windows_Bastide_ChrXtest.csv -cm /home/tiago/Documents/PoolLab/coldSIBSAM/Bastide_Window_cMorg.txt -o cold_args
# tiaaagosr@gmail.com, tsr@usp.br

import argparse
from vcfSite_class import *
from matedPair_class import *
from SAMread_class import *
from parentalDict_class import *
from matedPair_class import *
from mate_read_list_class import *
from window_class import *
from window_pair_class import *
from cMorg_window_info_class import *

from window_merge_function import *

min_read = 200
indel_remove_bp = 3

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chrm", help="Chromosome arm being studied.", required=True)
parser.add_argument("-p1", "--parent1", help="vcf file from the parent with the target phenotype (full path) (eg. cold parent in cold experiment).", required=True)
parser.add_argument("-p2", "--parent2", help="vcf file from the parent with the non-target phenotype (full path).", required=True)
parser.add_argument("-i1", "--indel1", help="INDEL vcf file from the parent with the target phenotype (full path) (eg. cold parent in cold experiment).", required=True)
parser.add_argument("-i2", "--indel2", help="INDEL vcf file from the parent with the non-target phenotype (full path).", required=True)
parser.add_argument("-f1", "--offspring1", help="samfile prefix (before chrm) from the offspring pool with the target phenotype (full path)", required=True)
parser.add_argument("-f2", "--offspring2", help="samfile prefix (before chrm) from the offspring pool with the non-target phenotype (full path)", required=True)
parser.add_argument("-w","--windows",help="genomic window tsv (start and end only).", required=True)
parser.add_argument("-cm","--cmorg",help="file with cMorg information for the windows (tab delimited, Chrm, Start, End, end-cMorg).", required=True)
parser.add_argument("-o","--output",help="output file name prefix (without file extension).", required=True)

args = parser.parse_args()

output_prefix = args.output
chrmCode = args.chrm
cMorg_filepath = args.cmorg
window_file = args.windows
parent1_filepath = args.parent1 # targetted phenotype parent (ie. cold ancestry parent in cold exp, highland/dark parent in pigm exp)
parent2_filepath = args.parent2
parent1_indel_filepath = args.indel1
parent2_indel_filepath = args.indel2
off1_sam_filepath = args.offspring1 # targetted phenotype pool (ie. cold pool in the cold experiment, dark pool in the pigm exp)
off2_sam_filepath = args.offspring2

if chrmCode == "4":
    chrmID = "X"
elif chrmCode == "3":
    chrmID = "2L"
elif chrmCode == "7":
    chrmID = "2R"
elif chrmCode == "5":
    chrmID = "3L"
elif chrmCode == "8":
    chrmID = "3R"
else:
    print("Chrm code and ID not recognized for the 5 chrm arms of Drosophila melanogaster.\n Edit the script accordingly.\n Quittind.\n")
    quit()

outlog = open(output_prefix+"_"+chrmID+".log","w")

## Using vcfSite class and parentalDict class
print("Calculating parental SNPs without indels.")
parent_indel_rm_pos = []
with open(parent1_indel_filepath) as p_indel:
    for r in p_indel:
        indel = vcfSite(r)
        indel_pos = int(indel.pos)
        indel_len = len(indel.ref)
        for rm_p in range(indel_pos-indel_remove_bp, indel_pos+indel_remove_bp+indel_len):
            parent_indel_rm_pos.append(str(rm_p))

with open(parent2_indel_filepath) as p_indel:
    for r in p_indel:
        indel = vcfSite(r)
        indel_pos = int(indel.pos)
        indel_len = len(indel.ref)
        for rm_p in range(indel_pos-indel_remove_bp, indel_pos+indel_remove_bp+indel_len):
            parent_indel_rm_pos.append(str(rm_p))
parent_indel_rm_pos = list(set(parent_indel_rm_pos))

parental_reads = parentalDict(parent1_filepath, parent2_filepath)
noIndel_parentalDict = {k: v for k, v in parental_reads.parental_dict.items() if k not in parent_indel_rm_pos}
outlog.write("Number of SNPs distinguinshing the parental genomes.\nWithout INDEL filtering: {}\nWith INDEL filtering: {}\n\n".format(len(parental_reads.parental_dict.keys()),len(noIndel_parentalDict.keys())))
print("Calculating informative reads for Pool 1.")
## mated_dictionary is a list of 4 values: 0: the dict itself, 1: total mated read pairs, 2: informative mated read pairs. 3: average number of SNPs per total read pairs. 1, 2, and 3 are useful to check data quality.
mated_reads, total_r, info_r, avg_snps = mate_read_list(off1_sam_filepath, noIndel_parentalDict).sorted_list()

print("Printing Pool1 log information.")
perc_read = round(info_r/total_r*100,2)
outlog.write("Pool 1!\nAverage number of SNPs per mated reads pair: {}.\nPercentage of informative reads = {}%. \nTotal mated read pairs: {}. Total informative read pairs: {}\n\n".format(avg_snps,perc_read, total_r, info_r)) 

print("Calculating informative reads per window for pool 1.")
## window
windowList1 = []
window_track = 0
with open(window_file) as wf:
    for r in wf:
        window_track += 1
        print(f"\rWindow: {window_track}", end="", flush=True)
        w = window(r)
        pool_ancestry1, pool_ancestry2, mated_reads = w.pool_ancestry(mated_reads)
        windowList1.append([w.start, w.end] + [pool_ancestry1, pool_ancestry2])

# delete previous Pool dict and make a new one for the 2nd pool
del mated_reads
print("\nCalculating informative reads for Pool 2.")
## mated_dictionary is a list of 4 values: 0: the dict itself, 1: total mated read pairs, 2: informative mated read pairs. 3: average number of SNPs per total read pairs. 1, 2, and 3 are useful to check data quality.
mated_reads, total_r, info_r, avg_snps = mate_read_list(off2_sam_filepath, noIndel_parentalDict).sorted_list()
print("Printing Pool2 log information.")
perc_read = round(info_r/total_r*100,2)
outlog.write("Pool 2!\nAverage number of SNPs per mated reads pair: {}.\nPercentage of informative reads = {}%. \nTotal mated read pairs: {}. Total informative read pairs: {}\n\n".format(avg_snps,perc_read, total_r, info_r)) 

print("Calculating informative reads per window for pool 2.")
windowList2 = []
window_track = 0
with open(window_file) as wf:
    for r in wf:
        w = window(r)
        window_track += 1
        print(f"\rWindow: {window_track}", end="", flush=True)
        pool_ancestry1, pool_ancestry2, mated_reads = w.pool_ancestry(mated_reads)
        windowList2.append([w.start, w.end] + [pool_ancestry1, pool_ancestry2])

print("\nCombining window information from both pools.")
window_full = []
if len(windowList1) == len(windowList2):
    for i in range(0,len(windowList1)):
        window_full.append(windowList1[i]+windowList2[i][-2:])
else:
    print("The two pools resulted in lists of different sizes.\nQuitting the script!!\n")
    quit()

print("Merging windows based on min_read_threshold: "+str(min_read)+".")
outlog.write("Minimum read threshold for merging windows was: "+str(min_read)+".\n")
merged1_window = []
for w in window_merge(window_full, min_read):
    if w[-1] != "NA":
        merged1_window.append(w)

# need to run twice in case the last window was preceeded by "NA" and did not pass min threshold itself
merged2_window = []
for w in window_merge(merged1_window, min_read):
    if w[-1] != "NA":
        merged2_window.append(w)

print("Calculating and printing ancestry difference window file for SIBSAM.")
ad_windows = []
for w in merged2_window:
    wpair = window_pair(w[0], w[1], w[2], w[3], w[4], w[5])
    p1f = wpair.pool1full
    p2f = wpair.pool2full
    ad = wpair.ancestry_difference
    ad_row = [chrmID]+[int(w[0]),int(w[1])]+[p1f, p2f, ad]
    ad_windows.append(ad_row)

# adding cMorg position to match the last base of the window
output = open(output_prefix+"_windows.tsv","w")
for w in ad_windows:
    wStart = int(w[1]) - 1 # removing the 1-base position that was used throughout the script because SIBSAM is 0-based
    wEnd = int(w[2]) - 1
    with open(cMorg_filepath) as cM_f:
        next(cM_f)
        for r in cM_f:
            cm = cMorg_window_info(r)
            if cm.chrm == w[0] and cm.end == wEnd:
                cMorg = cm.cMorg
                break
        new_w = [w[0]] + [wStart, wEnd, cMorg] + w[3:]
        outting_row = list(map(str,new_w))
        output.write("\t".join(outting_row)+"\n")
output.close()

outlog.write("windows used: "+window_file+"\n")
outlog.write("parent1 used: "+parent1_filepath+"\n")
outlog.write("parent2 used: "+parent2_filepath+"\n")
outlog.write("f1 pool1 used: "+off1_sam_filepath+"\n")
outlog.write("f1 pool2 used: "+off2_sam_filepath+"\n")

outlog.close()

