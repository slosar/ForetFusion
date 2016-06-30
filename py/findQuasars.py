#!/usr/bin/env python
##
## A very simple example file for Jose to follow
##
import numpy as np
import pyfits
print "Loading..."
fi=pyfits.open ('data//spAll-v5_10_0.fits')
da=fi[1].data
print "done"

# find matching masks for BOSS
bt1=[10,11,12,13,14,15,16,17,18,19,40,41,42,43,44] ## these are the qso target mask bits 
w_qso = np.zeros(len(da),dtype=bool)
for b in bt1:
    w_qso = w_qso | ((da.BOSS_TARGET1 & 2**b) > 0)
print "Found ",w_qso.sum(), "BOSS objects."

# find matching masks for eBOSS
ebt1=[10,11,12,13,14,15,16,17,18] ## these are the qso target mask bits 
ew_qso = np.zeros(len(da),dtype=bool)
for b in bt1:
    ew_qso = ew_qso | ((da.EBOSS_TARGET1 & 2**b) > 0)
    ew_qso = ew_qso | ((da.EBOSS_TARGET2 & 2**b) > 0)

print "Found ",ew_qso.sum(), "eBOSS objects."

w_tot=w_qso|ew_qso
print "Found", len(frozenset(da[w_qso].THING_ID)), "unique BOSS objects."
print "Found", len(frozenset(da[ew_qso].THING_ID)), "unique eBOSS objects."
print "Found", len(frozenset(da[w_tot].THING_ID)), "unique total objects."

#thing id dictionary
print "Copying..."
for tid,mjd,plate,fiber,ra,dec in zip(da[w_tot].THING_ID, da[w_tot].MJD, da[w_tot].PLATE,
                                      da[w_tot].FIBERID, da[w_tot].RA, da[w_tot].DEC):
    li.append((tid,(mjd,plate,fiber,ra,dec)))
print "Making dictionary..."
tid_dic=dict(li)

print "Dictionary length:", len(ti.keys())






