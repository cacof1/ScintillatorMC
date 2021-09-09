from ROOT import *
import sys
from root_numpy import hist2array
import numpy as np
from scipy.io import savemat
import matplotlib.pyplot as plt
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
"""
###

###
A = np.zeros((NX, NY),dtype=np.float32)
for i in range(0,NPB):
    Edep       = f.Get("YZProj_Q/YZProj_Q_"+str(i))
    Edep       = hist2array(Edep)
    A += Edep
A = np.flipud(A)
plt.imshow(A,cmap='gray_r')
plt.show()
"""
hist2 = f.Get("Front")
hist2 = hist2array(hist2)
hist2 = np.rot90(hist2)
plt.imshow(hist2,cmap='gray')
plt.show()
