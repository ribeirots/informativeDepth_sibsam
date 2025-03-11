# OOP implementation of the informative read depth call step for SIBSAM
# Class: vcfSites
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

import re

class matedPair:
    def __init__(self, votes):
        self.votes = votes

    @property
    def maxSNP(self):
        maxpos = self.votes[0][0]
        for v in self.votes:
            if v[0] > maxpos:
                maxpos = v[0]
        return maxpos

    @property
    def minSNP(self):
        minpos = self.votes[0][0]
        for v in self.votes:
            if v[0] < minpos:
                minpos = v[0]
        return minpos

    @property
    def ancestry(self):
        # votes based on genotype fully supporting one ancestry
        ANC1 = 0
        ANC2 = 0
        # votes based on genotype partially supporting one ancestry
        anc1 = 0
        anc2 = 0
        for v in self.votes:
            if v[1] == 1:
                ANC1 += v[1]
            elif v[1] == 0.5:
                anc1 += v[1]
            if v[2] == 1:
                ANC2 += v[2]
            elif v[2] == 0.5:
                anc2 += v[2]
        # lower case = half vote
        if ANC1 > ANC2:
            return "P1"
        elif ANC2 > ANC1:
            return "P2"
        elif anc1 > anc2:
            return "p1"
        elif anc2 > anc1:
            return "p2"
        else:
            return "tied"
