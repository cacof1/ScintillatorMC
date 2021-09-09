from ROOT import *
import numpy as np
import sys
from root_numpy import hist2array
import matplotlib.pyplot as plt
f=TFile(sys.argv[1])
print(f.ls())

hist = f.Get("Front")
print(dir(hist))
hist = hist.ProjectionXY("")
hist = hist2array(hist)
hist = np.rot90(hist)
plt.imshow(hist,cmap='gray')
plt.axis('off')
plt.savefig(sys.argv[1][:-5]+".png")
#np.savetxt("YZProj_Tot_Q.txt",hist)
