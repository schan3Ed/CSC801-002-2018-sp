import random
import time
import numpy as np
import sys
import math

BKV = 3.14159265359

def timedfunction(f, *a, **b):
    start = time.time()
    a = f(*a, **b)
    end = time.time()
    return end - start, a


def singleDrop(d=1.0, L=1.0):
    "dropping a single needle"
    y = np.random.uniform(0, d)
    angle = np.random.uniform(0, math.pi)
    height = L/2 * np.sin(angle)
  #  print y - height
    if (y + height) >= d or (y - height) <= 0:
        return True
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
            hit += 1.0
        cnt += 1.0
        if hit:
            pi = 2 * cnt / hit
        #print pi
        if pi >= BKV_lower and pi <= BKV_upper:
            return pi, cnt, False
    return pi, cnt, True

#@timedfunction
def run(probLmt=10000, tol=0.005, experimentCnt=5, seed=None):
    "main method for parallel line"
    entry = []
    seed = seed or np.random.randint(low=0, high=9999)
    np.random.seed(seed)
    for i in range(experimentCnt):
        isCensored = False
        t, result = timedfunction(singleExperiment, probLmt, tol, seed)
        pi, cnt, isCensored = result
        entry.append({
            "Pi": round(pi, 7), 
            "Count": cnt, 
            "IsCensored":isCensored, 
            "Seed":seed, 
            "Error": round(pi - BKV, 7)
            })
        seed = np.random.randint(low=0, high=9999)
        np.random.seed(seed)
    return entry

if __name__ == "__main__":
    p = run()
    for i in p:
        print i
  #  print("Program running")

