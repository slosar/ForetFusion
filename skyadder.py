#!/usr/bin/env python
from glob import glob
import numpy as np
import pyfits
import matplotlib.pyplot as plt
dirname="data/spectra_sky"
skip=10000
print "Preparing file list..."
filelist=glob(dirname+"/*/*.fits")
print "Done, we have ",len(filelist),"files."
nsp=0
loglam0=3.5 #3.45
loglamend=4.02 #4.1
loglamstep=1e-4
Np=(loglamend-loglam0)/loglamstep

meansky=np.zeros(Np)
meanw=np.zeros(Np)

Nc=3
Nit=5
contr=[np.zeros(Np) for i in range(Nc)]
loglam=loglam0+loglamstep*np.arange(Np)

for iter in range(Nit):
    for varcount in range(-1,Nc):
        nsp=0
        current=np.zeros(Np)
        currentw=np.zeros(Np)
        for ifile,filename in enumerate(filelist[::skip]):
            da=pyfits.open(filename)
            for ext in da[4:]:
                cloglam=ext.data["loglam"]
                csky=ext.data["sky"]+ext.data["flux"]
                civar=ext.data["ivar"]
                ## values of variables
                var=[np.cos(ext.header["ALT"]/180.*np.pi),
                            np.sin(ext.header["AZ"]/180.*np.pi),
                            np.cos(ext.header["AZ"]/180.*np.pi)]
                assert(len(var)==Nc)
                if (varcount==-1):
                    unwanted=np.zeros(Np)
                    for i in range(Nc):
                        unwanted+=contr[i]*var[i]
                else:
                    unwanted=meansky
                    for i in range(Nc):
                        if (i!=varcount):
                            unwanted+=contr[i]*var[i]

                ndx=np.array([int(v) for v in ((cloglam-loglam0)/loglamstep+0.5)])
		tmp = 1 if varcount == -1 else var[varcount]
		current[ndx]+=civar*tmp*(csky-unwanted[ndx])
                currentw[ndx]+=civar*tmp**2
                nsp+=1
        print current		
        current/=currentw
	if False:
	   plt.plot(current)
	   plt.ylim([0, 50])
	   plt.show()
	
        meanw=currentw
        if (varcount==-1):
            meansky=current
        else:
            contr[varcount]=current
	if False:
           plt.plot(current)
           plt.ylim([-50, 50])
           plt.show()
        print "Finished iteration/varcount",iter,varcount

    f=open("meansky%i.txt"%(iter),"w")
    for l,w,v in zip(loglam,meanw,meansky):
        if w>0:
            f.write("%g %g %g \n"%(l,w,v))
    f.close()
    for i in range(Nc):
        f=open("componet%i_%i.txt"%(i,iter),"w")
        for l,w,v in zip(loglam,meanw,contr[i]):
            if w>0:
                f.write("%g %g %g \n"%(l,w,v))
        f.close()
        
    

