import sys
import numpy as np
import scipy
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from root_numpy import hist2array
from ROOT import *
def errorFunction(x,a,b,c,d):
    return a*scipy.special.erf(b*(x-c))+d

def GetTheEdge(img,rot,xstart,xend,ystart,yend): #Function to compile the slainted edge into the ESF
    angle      = rot*np.pi/180.
    esf_res    = int(round(1./np.tan(angle)))
    print("Supersampled resolution:", esf_res, "New resolution", esf_res*(xend-xstart))
    print("Angle: ",angle)
    Xdata_res  = xend - xstart
    ESF      = np.zeros((Xdata_res*esf_res))
    Count    = np.zeros((Xdata_res*esf_res))
    xdata    = np.arange(0,Xdata_res,1)
    print(xdata)
    for id_ymin in range(ystart,yend,esf_res):
        ##Number of esf_ref we have done
        back_shift = int((yend-id_ymin)/esf_res)
        #back_shift = int((yend-id_ymin)*np.tan(angle))
        print(back_shift)
        for n in range(0, esf_res):
            if(n+id_ymin>=yend): break
            data = img[n+id_ymin, xstart:xend][::-1] #extract data from radiography
            #plt.plot(data)
            for idx,value in enumerate(data): ## Place at the appropriate value
                idbin =  esf_res*idx+ n - back_shift#*esf_res
                if(idbin>=Xdata_res*esf_res or idbin<0): continue
                ESF[idbin]+=value
                Count[idbin]+=1
    ESF[Count>0]/=Count[Count>0]
    return ESF

def LSF(edgespreadfunction): return np.diff(edgespreadfunction)

def MTF(lsf,dx): ## Function to take the discrete fourier transform of the LSF
    MTF_v = np.abs(np.fft.rfft(lsf))
    MTF_v /=MTF_v[0]
    xMTF_v = np.fft.rfftfreq(len(lsf),dx)
    return MTF_v,xMTF_v


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
tree = f.Get("Header")
data = np.zeros((NX,NY))
for event in tree: NPB = event.NPB
for i in range(0,NPB*NPB):
    hist = f.Get("YZProj_Q/YZProj_Q_"+str(i))
    hist = hist2array(hist)
    data +=hist

npixX = data.shape[1]
npixY = data.shape[0]

##Normalize the data
data  -= np.min(data)
data  /= np.max(data)
#data  = (data - np.mean(data))/np.std(data)

# pixel size
print("Image pixels in X/Y:", npixX, npixY)
pixSize = FOV/npixX

data    = np.fliplr(data)
xstart,xend, ystart, yend = 90, 120, 125,168
angle = 2.5

print("ROI start [xstart:xend,ystart:yend]",xstart,xend,ystart,yend)
###############################################################################
## Print full image slice, including ROI drawn in
#Draw image

plt.imshow(data,origin='lower', cmap=cm.gray)#, zorder=0, extent=[x[0], x[-1], y[0], y[-1]], vmin=0.1, vmax =2.1) #0.1, 2.1
#plt.plot([xstart,xend],[ystart,ystart],'r--')
#plt.plot([xstart,xend],[yend,yend],'r--')
#plt.plot([xstart,xstart],[ystart,yend],'r--')
#plt.plot([xend,xend],[ystart,yend],'r--')
plt.show()

###############################################################################
##Start with MTF calculation
##Get the Edge

ESF = GetTheEdge(data,angle,xstart,xend,ystart,yend) ###FIXME: This angle might not be the correct one!
Xdata_res = xend - xstart
XESF = np.linspace(-1.,Xdata_res-1, len(ESF),endpoint=True)
popt,pcov = opt.curve_fit(errorFunction, XESF, ESF)
ESFfit = errorFunction(XESF, *popt)
plt.plot(XESF, ESF, label="Meas. ESF")
plt.plot(XESF, ESFfit, label="Error-function fit")
plt.legend()
plt.show()

lsf = LSF(ESFfit)
mtf,xmtf = MTF(lsf,pixSize*np.tan(np.deg2rad(angle))) ##FIXME: Correct angle?

plt.plot(xmtf,mtf,label="MTF for cube at "+str(angle)+"deg.")
plt.legend()
plt.show()

print('MTF_10%=',xmtf[np.where(mtf<=0.1)][0],'lp/mm')
