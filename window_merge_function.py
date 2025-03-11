# OOP implementation of the informativssistire read depth call step for SIBSAM
# Designed to run 1 chrm arm at a time!
# tsr@usp.br

def window_merge(windowList, min_read_count):
    for i in range(0,len(windowList)):
        merge_by_pool1 = 0
        merge_forward = 0

        pool1reads = windowList[i][-4]+windowList[i][-3]
        pool2reads = windowList[i][-2]+windowList[i][-1]
        # if True: merge is necessary
        if pool1reads < min_read_count or pool2reads < min_read_count:
            # will perform merge based on which pool has fewer reads
            if pool1reads < pool2reads:
                merge_by_pool1 = 1
            # not the first or last window
            if i > 0 and i < len(windowList)-1:
                # if True, the previous window was already merged to a different one
                if windowList[i-1][-4] == "NA":
                    merge_forward = 1
                # if the previous window is valid, checks (based on smaller pool) which window (previous or next) has fewer reads and will receive these reads
                elif merge_by_pool1:
                    if windowList[i-1][-4] + windowList[i-1][-3] > windowList[i+1][-4] + windowList[i+1][-3]:
                        merge_forward = 1
                # merge by pool2
                else:
                    if windowList[i-1][-2] + windowList[i-1][-1] > windowList[i+1][-2] + windowList[i+1][-1]:
                        merge_forward = 1
            elif i == 0:
                merge_forward = 1
            
            if merge_forward:
                windowList[i+1][0] = windowList[i][0] # change start position of the next window
                windowList[i+1][-4] += windowList[i][-4]
                windowList[i+1][-3] += windowList[i][-3]
                windowList[i+1][-2] += windowList[i][-2]
                windowList[i+1][-1] += windowList[i][-1]
                windowList[i][-4:] = ["NA", "NA", "NA","NA"]
            # It is possible that in the last window it will try to merge back into an invalid window. This if prevents that, but requires that this loop needs to be ran 2x - in the raw window List and in the 1st merged window list to ensure that the last window also passes the filtering requirement. Otherwise, the last window could remain below min threshold
            elif windowList[i-1][-4] != "NA":
                windowList[i-1][1] = windowList[i][1] # change end position of the previous window
                windowList[i-1][-4] += windowList[i][-4]
                windowList[i-1][-3] += windowList[i][-3]
                windowList[i-1][-2] += windowList[i][-2]
                windowList[i-1][-1] += windowList[i][-1]
                windowList[i][-4:] = ["NA", "NA", "NA","NA"]
    return windowList

