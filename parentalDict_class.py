# OOP implementation of the informativssistire read depth call step for SIBSAM
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

from vcfSite_class import *

## Using vcfSite class
class parentalDict:
    def __init__(self, filepath1, filepath2):
        self.filepath1 = filepath1
        self.filepath2 = filepath2

    @property
    def parental_dict(self):
        with open(self.filepath1) as p1, open(self.filepath2) as p2:
            chrmDict = {}
            for line1, line2 in zip(p1, p2):
                # check if both files have a header of the same length, aborts if they don't
                if (line1[0] == "#" and line2[0] != "#") or (line1[0] != "#" and line2[0] == "#"):
                    print("Error 1: the parental vcf files have headers of different sizes!\nQuitting the script!!!")
                    quit()
                elif line1[0] != "#":
                    a = vcfSite(line1)
                    b = vcfSite(line2)
                    if a.isvalidsnp(b):
                        site_pos_allele = a.get_pos_alleles(b)
                        chrmDict[site_pos_allele[0]] = [site_pos_allele[1], site_pos_allele[2]]
        return chrmDict

    def __repr__(self):
        pass
    
    def __str__(self):
        pass

