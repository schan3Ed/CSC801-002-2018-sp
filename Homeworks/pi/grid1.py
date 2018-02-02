
BKV = 3.14

class grid1():
    def singleDrop(self, lw, nl):
        "dropping a single needle"
        return False
        
    def singleExperiment(self, lw, nl, pl, err):
        "running a single experiment"
        return

    def __main__(self, lineWidth=10, needleLength=20, probLmt=100, 
        errInterval=0.05, experimentCnt=50):
        "main method for parallel line"

        for i in range(experimentCnt):
            singleExperiment(lineWidth, needleLength, probLmt, errInterval)
        return

