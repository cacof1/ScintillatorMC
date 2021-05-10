import numpy as np
from ROOT import *
import matplotlib.pyplot as plt
import sys
from root_numpy import hist2array
from scipy.io import savemat

A = TFile(sys.argv[1])
print(A.ls())
hist = hist2array(A.Get("L_Tot"))
dc   = {"hist":hist}
savemat("matlab_matrix.mat", dc)
