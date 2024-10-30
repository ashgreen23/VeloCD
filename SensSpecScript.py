import sys
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects

utils = importr('utils')

path=str(sys.argv[10])
def RunSensSpecificity(ResultFileMyMetadata, GeneSignature, MyMetadata, MyMethod, Output, Class1, Class2, InterName, PredGroup, thresholdProb):
    r=ro.r
    r.source(path+"SensitivitySpecificity.R")
    p=r.SensitivitySpecificity(ResultFileMyMetadata, GeneSignature, MyMetadata, MyMethod, Output, Class1, Class2, InterName, PredGroup, thresholdProb)
    return p

a=RunSensSpecificity(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]), str(sys.argv[7]), str(sys.argv[8]), str(sys.argv[9]), float(sys.argv[11]))
