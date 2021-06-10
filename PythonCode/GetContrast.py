import sys
import numpy as np
import scipy
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from root_numpy import hist2array
from ROOT import *
def create_circular_mask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

try:
    fileName = sys.argv[1]
    FOV = float(sys.argv[2])
except:
    print( 'Error. Use', sys.argv[0], ' [txtFileName] [FOV/mm] ')
    exit

###############################################################################
## Get the image data
f=TFile(fileName)
print(f.ls())

f=TFile(sys.argv[1])
hist = hist2array(f.Get("YZProj_Q/YZProj_Q_0"))
NX, NY = hist.shape
NPB = 0

Xmin, Xmax = -FOV/2, FOV/2
tree = f.Get("Header")
data = np.zeros((NX,NY))
for event in tree: NPB = event.NPB
for i in range(0,NPB):
    hist = f.Get("YZProj_Q/YZProj_Q_"+str(i))
    hist = hist2array(hist)
    temp = hist/np.max(hist)
    data +=hist    

## Forms a block of 6x5
PosX      =np.array([-75,-45,-15,15,45,75])  # mm - 6
PosY      =np.array([-60,-30,0,30,60])  #mm - 5
Diameter  =np.array([0.5/2, 2./2, 4./2, 7./2,   10./2, 15./2])  #mm -6
Depth     =np.array([0.5/2, 1./2, 2./2, 3.25/2, 4.5/2])  #mm -5

## Idem in the image space
PosX      = (PosX - Xmin)*(NX/FOV)# pix - 6
PosY      = (PosY - Xmin)*(NX/FOV)# pix - 5
Diameter  = Diameter*(NX/FOV) #pix -6
Depth     = Depth*(NX/FOV)  #pix -5
out       = []
mask_tot = np.zeros((150,150))
for idx,x in enumerate(PosX): ## Loop over all holes
    for idy,y in enumerate(PosY):

        mask = create_circular_mask(150,150,center = [x, y],radius = Diameter[idx])
        mask_tot += mask
        mean, sigma = np.average(data[mask]), np.std(data[mask])
        out.append([Diameter[idx],Depth[idy], mean,sigma])


## Now for the background
mask = create_circular_mask(150,150,center = [40,40],radius = 5)
mask_tot += mask
plt.imshow(data)
plt.imshow(mask_tot,alpha=0.5)
plt.show()
mean, sigma = np.average(data[mask]), np.std(data[mask])
out.append([5,0, mean,sigma])
np.savetxt("out.txt",out)


## For when it is converted to WET
#RSP_al = 2.16282
#RSP_w  = 1.0
#WETDiff = Depth*(RSP_al-RSP_w)
