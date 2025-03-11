# OOP implementation of the informative read depth call step for SIBSAM
# Class: vcfSites
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

import re

class SAMread:
    n_indel_rm = 3
    
    def __init__(self, samstr):
        self.samlist = samstr[:-1].split("\t")
        self.name = self.samlist[0]
        self.flag = self.samlist[1]
        self.chrm = self.samlist[2]
        self.pos = self.samlist[3]
        self.qual = self.samlist[4]
        self.cigar = self.samlist[5]
        self.matestart = self.samlist[7]
        self.seq = self.samlist[9]


    @property
    def cigar_split(self):
        cigar_list = list(re.findall(r'(\d+)([a-zA-Z])', self.cigar))
        if "I" not in self.cigar or len(cigar_list) == 1:
            return [[cig_code, int(cig_count)] for cig_count, cig_code in cigar_list]
        else:
            indel_cigar_list = []
            for i in range(0, len(cigar_list)):
                cigar_list[i] = list(cigar_list[i])
                cigar_list[i][0] = int(cigar_list[i][0])
            # if there is at least one indel, loop through cigar string list
            for i in range(0, len(cigar_list)):
                # if not in the last value
                if i < len(cigar_list)-1:
                    if cigar_list[i][1] == "M":
                        if cigar_list[i+1][1] == "I" or cigar_list[i+1][1] == "D":
                            if cigar_list[i][0] - self.n_indel_rm >= 0:
                                new_cigar = [self.n_indel_rm, "R"]
                                cigar_list[i][0] -= self.n_indel_rm
                            else:
                                new_cigar = [cigar_list[i][0], "R"]
                                cigar_list[i][0] = 0
                            indel_cigar_list.append(cigar_list[i])
                            indel_cigar_list.append(new_cigar)
                        else:
                            indel_cigar_list.append(cigar_list[i])
                    elif cigar_list[i][1] == "I" or cigar_list[i][1] == "D":
                        indel_cigar_list.append(cigar_list[i])
                        if cigar_list[i+1][1] == "M":
                            if cigar_list[i+1][0] - self.n_indel_rm >= 0:
                                new_cigar = [self.n_indel_rm, "R"]
                                cigar_list[i+1][0] -= self.n_indel_rm
                            else:
                                new_cigar = [cigar_list[i+1][0], "R"]
                                cigar_list[i+1][0] = 0
                            indel_cigar_list.append(new_cigar)
                    else:
                        indel_cigar_list.append(cigar_list[i])
                else:
                    indel_cigar_list.append(cigar_list[i])
            return [[cig_code, int(cig_count)] for cig_count, cig_code in indel_cigar_list]

    def query_ref_pos(self):
        query_i = 0
        ref_i = int(self.pos)
        read_bases = []
        for cigar_step in self.cigar_split:
            if cigar_step[0] == "M" or cigar_step[0] == r"=" or cigar_step[0] == "X": # consumes query and reference
                for i in range(0,cigar_step[1]):
                    read_bases.append([str(ref_i), self.seq[query_i]])
                    query_i += 1
                    ref_i += 1
            elif cigar_step[0] == "I" or cigar_step[0] == "S": # consumes query but not reference
                for i in range(0,cigar_step[1]):
                    query_i += 1
            elif cigar_step[0] == "D" or cigar_step[0] == "N": # consumes reference but not query
                for i in range(0,cigar_step[1]):
                    ref_i += 1
            # I created the "R" cigar code to include regions to "Remove" with filters
            # The filter being applied is a indel filter - removing 3bp around indel.
            # It will walk on both the query and reference as these are applied to mapped regions, but ignores them
            elif cigar_step[0] == "R": # consumes reference and query WITHOUT reading bases
                for i in range(0,cigar_step[1]):
                    query_i += 1
                    ref_i += 1
            elif cigar_step[0] != "H" and cigar_step[0] != "P": # H and P consumes nothing, if it is another code: shows error
                print("class SAMread Error 1: unknown cigar code detected.\nQuitting script.\n")
        return read_bases

    # parentalDict KEYS are positions of SNPs, and and VALUE is a list of two lists with the genotypes of each parent
    # return output: list of vote per position, [position, voteParent1, voteParent2]
    def ancestry_votes(self, parentalDict):
        readquery = self.query_ref_pos()
        ancVotes = []
        for base in readquery:
            if str(base[0]) in parentalDict.keys():
                parent1, parent2 = parentalDict[str(base[0])]
                if base[1] in parent1 and base[1] not in parent2:
                    ancVotes.append([base[0],1,0])
                elif base[1] in parent1: # will also be in parent2
                    if parent1.count(base[1]) == 2 and parent2.count(base[1]) == 1:
                        ancVotes.append([base[0],0.5,0])
                    elif parent1.count(base[1]) == 1 and parent2.count(base[1]) == 2:
                        ancVotes.append([base[0],0,0.5])
                    else:
                        print("class SAMread Error 2: Genotypes not expected by the script. See offspring/sam genotype below followed by parental genotypes from VCF and dictionary.\nQuitting script!!!\n")
                        print(base[1])
                        print(parent1)
                        print(parent2)
                        quit()
                elif base[1] in parent2: # will not be in parent1
                    ancVotes.append([base[0],0,1])
                else: # it is not present in any parent
                    pass
        return ancVotes

    def __repr__(self):
        return 'SAMread("{}")'.format(r"\t".join(self.samlist))

    def __str__(self):
        return "SAMread - Name:{}, Chrm:{}, Pos:{}.".format(self.name, self.chrm, self.pos)

