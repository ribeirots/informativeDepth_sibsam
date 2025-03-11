# OOP implementation of the informative read depth call step for SIBSAM
# Class: vcfSites
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

class vcfSite:
    qual_threshold = 20
    min_sample_depth = 1

    def __init__(self, vcfstr): # receives the entire row from vcf file as a string
        self.vcflist = vcfstr[:-1].split("\t")
        self.chrm = self.vcflist[0]
        self.pos = self.vcflist[1]
        self.ref = self.vcflist[3]
        self.alt = self.vcflist[4]
        self.qual = self.vcflist[5]
        self.vcfformat = self.vcflist[8]
        self.vcfvalue = self.vcflist[9]
    
    # Gets the index of the value containing the GT information (genotype)
    @property
    def GTindex(self):
        formatcontent = self.vcfformat.split(":")
        if "GT" in formatcontent:
            return formatcontent.index("GT")
        else:
            return None
 
    # Gets the index of the value containing the DP information (total depth)
    @property
    def DPindex(self):
        formatcontent = self.vcfformat.split(":")
        if "DP" in formatcontent:
            return formatcontent.index("DP")
        else:
            return None

    # Gets the the value containing the GT information (genotype)
    @property
    def GTvalue(self):
        thisGT = self.vcfvalue.split(":")[self.GTindex]
        if "/" in thisGT:
            return thisGT.split("/")
        elif "|" in thisGT:
            return thisGT.split("|")
        else:
            print("Error 2: non-diploid site.")
            quit()
            return None
    
    # Gets the the value containing the DP information (total depth)
    @property
    def DPvalue(self):
        return self.vcfvalue.split(":")[self.DPindex]
    
    # Checks if str length was 10, if it wasn't it is an indication that maybe the vcf was generated from multiple or merged bam files and will be shown as multisample. In this case, the logics for genotyping this specific file will not hold and will need class will need to be adjusted with code to either pick one of the samples or combine their information.
    def isonesample(self):
        if len(self.vcflist) == 10:
            return True
        else:
            return False
    
    # Checks if the quality of the site is higher than the set quality threshold
    def isqual(self):
        if self.qual != '.':
            if float(self.qual) >= self.qual_threshold:
                return True
            else:
                return False
        else:
            return False
    
    # Checks whether two sites are in the same position
    def issamepos(self, pair):
        if self.chrm == pair.chrm and self.pos == pair.pos:
            return True
        else:
            return False

    # Checks whether the site pass min depth threshold
    def ismindepth(self):
        if int(self.DPvalue) >= self.min_sample_depth:
            return True
        else:
            return False

    # Checks if site pair has exactly two alleles (if fewer it is not a SNP, it higher it is not biallelic)
    # Also excludes sites that have two alleles and are identical heterozygous for both self and pair
    def isbiallelicSNP(self, pair):
        if len(set(self.getalleles() + pair.getalleles())) == 2 and set(self.getalleles()) != set(pair.getalleles()):
            return True
        else:
            return False

    # Checks if the site has a diploid genotype
    def isdiploid(self):
        if len(self.GTvalue) == 2:
            return True
        else:
            return False

    # Get the nucleotidic base code for the genotype of the site
    def getalleles(self):
        if len(self.alt) == 1:
            if self.GTvalue == ["0","0"]:
                return [self.ref, self.ref]
            elif self.GTvalue == ["1","1"]:
                return [self.alt, self.alt]
            else:
                return [self.ref, self.alt]
        elif len(self.alt.split(",")) == 2:
            return [self.alt.split(",")[0], self.alt.split(",")[1]]
        else:
            print("Error 3. Investigate the following Alt genotype further and fix this code.\n")
            print(self.alt)
            quit()
            return None

    # Checks if the site pair pass all the criteria to be considered a SNP for the target chromosome arm
    def isvalidsnp(self, pair):
        if self.issamepos(pair):
            if self.isonesample() and pair.isonesample():
                if self.isdiploid() and pair.isdiploid():
                    if self.isqual() and pair.isqual():
                        if self.ismindepth() and pair.ismindepth():
                            if self.isbiallelicSNP(pair):
                                return True
                    return False
                else:
                    print("Error 4. At least one of the VCF files is not diploid and the script needs to be adjusted to handle that.\nQuitting script!!!")
                    print(self.vcflist)
                    print(pair.vcflist)
                    quit()
                    return None
            else:
                print("Error 5. At least one of the VCF files contain multiple samples and the script needs to be adjusted to handle that.\nQuitting script!!!")
                print(self.vcflist)
                print(pair.vcflist)
                quit()
                return None
        else:
            print("Error 6. The files do not match number of rows and positions.\nQuitting script!!!")
            print(self.vcflist)
            print(pair.vcflist)
            quit()
            return None
    def printalleles(self, pair):
        return [self.getalleles(), pair.getalleles()]
    
    def get_pos_alleles(self, pair):
        return [str(self.pos), self.getalleles(), pair.getalleles()]

    def __repr__(self):
        return 'vcfSite("{}")'.format(r"\t".join(self.vcflist))

    def __str__(self):
        return "VCF site - Chrm:{}, Pos:{}, Ref:{}, Alt:{}, Genotype:{}".format(self.chrm, self.pos, self.ref, self.alt, self.GTvalue)



