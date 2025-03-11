# OOP implementation of the informativssistire read depth call step for SIBSAM
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

#from vcfSite_class import *
from SAMread_class import *
from matedPair_class import *
#from parentalDict_class import *


class mate_read_list:
    def __init__(self, sampath, parental_dict):
        self.sampath = sampath
        self.parental_dict = parental_dict

    @property
    def mated_dictionary(self):
        matedDict = {}
        with open(self.sampath) as samf:
            for r in samf:
#                if len(matedDict.keys()) == 30:
#                    break
                if r[0] != r"@":
                    current_read = SAMread(r)
                    if current_read.name in matedDict.keys():
                        if len(matedDict[current_read.name]) == 1:
                            matedDict[current_read.name].append(current_read.ancestry_votes(self.parental_dict))
                        else:
                            print("matedDict key has two reads already ans is trying to add a third. Supplementary and secondary aligments should've been removed prior to running this script, double check the filter steps and alter this script if needed.\nQuitting script!!!\n")
                            print(matedDict)
                            quit()
                    else:
                        matedDict[current_read.name] = [current_read.ancestry_votes(self.parental_dict)]

        unmated_keys = []
        uninformative_keys = []
        total_info_snps = 0
        for k in matedDict.keys():
            # del key if there is only a single read instead of mated pair
            if len(matedDict[k]) == 1:
                unmated_keys.append(k)
            # del key if the mated pair has no informative read
            elif len(matedDict[k][0]) + len(matedDict[k][1]) == 0:
                uninformative_keys.append(k)

        for k in unmated_keys:
            del matedDict[k]
        total_mated_reads = len(matedDict.keys())

        for k in uninformative_keys:
            del matedDict[k]

        tied_ancestry_keys = []
        for k in matedDict.keys():
            votes = []
            for v in matedDict[k][0] + matedDict[k][1]:
                if v not in votes:
                    votes.append(v)

            for i in range(0,len(votes)):
                votes[i][0] = int(votes[i][0])

            # search for SNPs in each each mated read had different votes for the same position
            snp_vote_check = []
            snp_del_vote = [] # list of SNPs with conflicting votes (likely sequecing error, each pair voted differently for the same position) that will be deleted
            for v in votes:
                if v[0] not in snp_vote_check:
                    snp_vote_check.append(v[0])
                else:
                    snp_del_vote.append(v[0])
            valid_votes = []
            for v in votes:
                if v[0] not in snp_del_vote:
                    valid_votes.append(v)
            read_pair = matedPair(valid_votes)
            total_info_snps += len(valid_votes)
            if read_pair.ancestry == "tied":
                tied_ancestry_keys.append(k)
            matedDict[k] = valid_votes

        for k in tied_ancestry_keys:
            del matedDict[k]
        informative_reads = len(matedDict.keys())

        return [matedDict, total_mated_reads, informative_reads, round(total_info_snps/total_mated_reads,2)]

    def sorted_list(self):
        mated_dict, total_reads, info_reads, av_snp = self.mated_dictionary
        mated_list = list(mated_dict.values())
        for i in range(0,len(mated_list)):
            mated_list[i] = sorted(mated_list[i], key=lambda x: x[0])
        mated_list = sorted(mated_list, key=lambda x: x[0][0])

        return [mated_list, total_reads, info_reads, av_snp]



