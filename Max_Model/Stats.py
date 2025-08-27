import numpy as np

# for ALL functions predicted and actual are sets
def calcFalsePositive(predicted, actual):
    fp = predicted - actual
    return len(fp)
def calcTruePositive(predicted, actual):
    tp = predicted & actual
    return len(tp)
def calcFalseNegative(predicted, actual):
    fn = actual - predicted
    return len(fn)
def calcTrueNegative(predicted, actual, numCells):
    tn = numCells - len(predicted|actual)
    return tn

# fp is false positive, tp is true positive, fn is false negative, fp is false positive. All are ints
def calcAccuracy(fp, fn, tp):
    acc = (tp)/(tp+fp+fn)
    return acc
def calcSpecificity(fp, tn):
    spec = tn/(tn+fp)
    return spec
def calcPrecision(fp, tp):
    prec = tp/(tp+fp)
    return prec
def calcRecall(fn, tp):
    recall = tp/(tp+fn)
    return recall

def calcKappa(fp, tp, fn):
    if fp != 0 and fn != 0:
        k = (2*((tp)-(fn*fp)))/((tp+fp)*(fp)*(tp+fn)*(fn))
    elif fp != 0 and fn == 0:
        k = (2*((tp)-(fn*fp)))/((tp+fp)*(fp)*(tp+fn))
    elif fn != 0 and fp == 0:
        k = (2*((tp)-(fn*fp)))/((tp+fp)*(tp+fn)*(fn))
    elif fp == 0 and fn == 0:
        k = (2*((tp)-(fn*fp)))/((tp+fp)*(tp+fn))
    return k

def calcQuantityDisagreement(predicted, actual):
    qd = np.abs(len(predicted) - len(actual))
    return qd
def calcAllocationDisagreement(fp, fn, qd):
    ad = fp + fn - qd
    return ad