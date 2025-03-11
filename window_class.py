# OOP implementation of the informative read depth call step for SIBSAM
# Class: window
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

from matedPair_class import * 

class window:
    def __init__(self, windowstr): # receives the entire row from vcf file as a string
        self.windowlist = windowstr[:-1].split("\t")
        self.start = int(self.windowlist[0])+1 # adding one because VCF and SAM are 1-based
        self.end = int(self.windowlist[1])+1

    def pool_ancestry(self, read_list):
        parent1 = 0
        parent2 = 0
        for i in range(0,len(read_list)):
            readpair = read_list[i]
            k_pair = matedPair(readpair)
            firstSNP = k_pair.minSNP
            lastSNP = k_pair.maxSNP
            # read after the window is over.
            # This only works because read list is ordered by firstSNP
            if firstSNP > self.end:
                break
            # read pair fully within the window
            elif firstSNP >= self.start and lastSNP <= self.end: 
                if k_pair.ancestry == "P1":
                    parent1 += 1
                    read_list[i].append("done")
                elif k_pair.ancestry == "P2":
                    parent2 += 1
                    read_list[i].append("done")
                elif k_pair.ancestry == "p1":
                    parent1 += 0.5
                    read_list[i].append("done")
                elif k_pair.ancestry == "p2":
                    parent2 += 0.5
                    read_list[i].append("done")
            # read part partially within
            elif firstSNP <= self.end and lastSNP >= self.start:
                window_votes = 0
                for v in k_pair.votes:
                    if int(v[0]) >= self.start and int(v[0]) <= self.end:
                        window_votes += 1
                # if most voting SNPs are in this window, or exact half the votes including the first SNP, count the read pair as belonging to this window
                prop_votes = window_votes / len(k_pair.votes)
                if (prop_votes > 0.5) or (prop_votes == 0.5 and firstSNP >= self.start and firstSNP <= self.end):
                    if k_pair.ancestry == "P1":
                        parent1 += 1
                        read_list[i].append("done")
                    elif k_pair.ancestry == "P2":
                        parent2 += 1
                        read_list[i].append("done")
                    elif k_pair.ancestry == "p1":
                        parent1 += 0.5
                        read_list[i].append("done")
                    elif k_pair.ancestry == "p2":
                        parent2 += 0.5
                        read_list[i].append("done")
        # remove reads that were already analyzed from read_list
        read_list = [x for x in read_list if x[-1] != "done"]
        return [parent1, parent2, read_list]

