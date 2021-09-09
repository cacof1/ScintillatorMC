from ROOT import *
import sys
from root_numpy import hist2array
import numpy as np
from scipy.io import savemat
f = TFile(sys.argv[1])
print(f.ls())



##
##
## First getting the info
hist = hist2array(f.Get("YZProj_Q/YZProj_Q_0"))
NX, NY = hist.shape
NPB = 0
tree = f.Get("Header")
for event in tree: NPB = event.NPB

###

###
A = np.zeros((NPB, NX, NY),dtype=np.float32)
for i in range(0,NPB):
    Edep       = f.Get("YZProj_Q/YZProj_Q_"+str(i))
    Edep       = hist2array(Edep)
    A[i] = Edep
savemat(sys.argv[1][:-5]+".mat",{"arr":A})
