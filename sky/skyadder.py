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
loglam0=3.45
loglamend=4.1
loglamstep=1e-4
Np=(loglamend-loglam0)/loglamstep

meansky=np.zeros(Np)

Nc=2
Nit=50
contr=[np.zeros(Np) for i in range(Nc)]
loglam=loglam0+loglamstep*np.arange(Np)

dotest=True
loglammean=loglam0+Np/2*loglamstep
tmean=100*(loglam-loglammean)**2
tcomp1=np.ones(Np)
tcomp2=loglam
tcomp3=np.sin(loglam)


for iter in range(Nit):
    for varcount in range(-1,Nc):
        nsp=0
        current=np.zeros(Np)
        currentw=np.zeros(Np)
        for ifile,filename in enumerate(filelist[::skip]):
            print filename
            da=pyfits.open(filename)
            for ext in da[4:]:
                cloglam=ext.data["loglam"]
                csky=ext.data["sky"]+ext.data["flux"]
                civar=ext.data["ivar"]
                ndx=np.array([int(v) for v in ((cloglam-loglam0)/loglamstep+0.5)])
                ## values of variables
                var=[np.cos(ext.header["ALT"]/180.*np.pi)-0.5,
                    np.sin(ext.header["AZ"]/180.*np.pi),
                    np.cos(ext.header["AZ"]/180.*np.pi)]
                var=var[:Nc]
                
                assert(len(var)==Nc)
                if dotest:
                    csky=tmean[ndx]+tcomp1[ndx]*var[0]+tcomp2[ndx]*var[1]#+tcomp3[ndx]*var[2]
                if (varcount==-1):
                    unwanted=np.zeros(Np)
                    for i in range(Nc):
                        unwanted+=contr[i]*var[i]
                else:
                    unwanted=meansky*1.0
                    for i in range(Nc):
                        if (i!=varcount):
                            unwanted+=contr[i]*var[i]

                tmp = 1 if varcount == -1 else var[varcount]
                current[ndx]+=civar*tmp*(csky-unwanted[ndx])
                currentw[ndx]+=civar*tmp**2
                nsp+=1
        print(current, currentw)
        current/=currentw
        current[np.where(current>10)]=10.0
        current[np.where(current<0)]=0.0
        
        if (varcount==-1):

            meansky=current*1.0

            f=open("meansky%i.txt"%(iter),"w")
            for l,w,v in zip(loglam,currentw,meansky):
                print l,w,v
                if w>0:
                    f.write("%g %g %g %g \n"%(l,w,v,100*(l-loglammean)**2))
            f.close()

        else:
            contr[varcount]=current*1.0
            f=open("componet%i_%i.txt"%(varcount,iter),"w")
            for l,w,v in zip(loglam,currentw,current):
                if w>0:
                    if (varcount==0):
                        tr=1
                    elif varcount==1:
                        tr=l
                    elif varcount==2:
                        tr=np.sin(l)
                    f.write("%g %g %g %g  \n"%(l,w,v,tr))
            f.close()

        if False:
            plt.plot(current)
            plt.ylim([0, 50])
            plt.show()
            plt.plot(current)
            plt.ylim([-50, 50])
            plt.show()
    print "Finished iteration/varcount",iter,varcount
        
    

