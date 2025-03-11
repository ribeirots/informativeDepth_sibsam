# OOP implementation of the informative read depth call step for SIBSAM
# Class: window
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

class window_pair():
    def __init__(self,start,end,pool1p1,pool1p2,pool2p1,pool2p2):
        self.start = int(start)
        self.end = int(end)
        self.pool1p1 = int(pool1p1)
        self.pool1p2 = int(pool1p2)
        self.pool2p1 = int(pool2p1)
        self.pool2p2 = int(pool2p2)
        self.pool1full = pool1p1 + pool1p2
        self.pool2full = pool2p1 + pool2p2
    
    @property
    def ancestry_difference(self):
        if self.pool1full > 0 and self.pool2full > 0:
            return (self.pool1p1/self.pool1full) - (self.pool2p1/self.pool2full)
        else:
            return None
