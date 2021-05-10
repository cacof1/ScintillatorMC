import numpy as np
from ROOT import *
import matplotlib.pyplot as plt
import sys
from root_numpy import hist2array
from scipy.io import savemat

data = np.zeros((2500,150,150))
A = TFile(sys.argv[1])
for i in range(2500):
    print(A.Get("YZProj_Q/YZProj_Q_"+str(i)))
    hist = A.Get("YZProj_Q/YZProj_Q_"+str(i)))
    hist = hist2array(hist)
    data[i] = hist

#hist = hist2array(A.Get("L_Tot"))
dc   = {"hist":data}
savemat("Slantde_edge_matlab_matrix.mat", dc)
