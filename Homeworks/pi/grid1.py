import random
BKV = 3.14
def singleDrop():
    "dropping a single needle"
    return False
    
def singleExperiment(pl, err, seed):
    "running a single experiment"
    cnt = 0
    hit = 0
    pi = 0
    BKV_upper = err + BKV
    BKV_lower = BKV - err
    for i in range(pl):
        if(singleDrop()):
            hit += 1
        cnt += 1
        pi = 2 * hit / cnt
        if pi > BKV_lower and pi < BKV_upper:
            return pi, cnt, False
    return -1, cnt, True

def run(probLmt=100, errInterval=0.05, experimentCnt=50, seed=None):
    "main method for parallel line"
    entry = []
    seed = seed or random.random() * 9999
    for i in range(experimentCnt):
        isCensored = False
        pi, cnt, isCensored = singleExperiment(probLmt, errInterval, seed)
        entry.append([pi, cnt, isCensored])
    print(entry)
    return

if __name__ == "__main__":
    run()

