# OOP implementation of the informative read depth call step for SIBSAM
# Class: window
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

class cMorg_window_info():
    def __init__(self,cMorgstr):
        self.cMorglist = cMorgstr[:-1].split("\t")
        self.chrm = self.cMorglist[0]
        self.start = self.cMorglist[1]
        self.end = int(self.cMorglist[2])
        self.cMorg = self.cMorglist[3]

