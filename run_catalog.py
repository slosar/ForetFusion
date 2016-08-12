
"""
Running the qso_catalog.py:

* First it goes through the file and filter all the qsos in the BOSS, EBOSS
that satisfy the bit condition and also
    'CLASS== "QSO" & (OBJTYPE=="QSO" | ''OBJTYPE=="NA".) & THING_ID != -1'
* Computes the healpix given the RA and DEC, and group all of them by healpix.
* Given a healpix, finds the objects with the same THING_ID
* Coadd them (average ivar*flux.) and compute the chisq
* Loop over:  If the chisq is more than 4, eliminate that spec and
coadd again and get new chisq
"""

from qso_catalog import *
from get_files import *

Pars = Ini_params()

# read the full SpAll file or the subset we're interested on
if False: df_fits = Ini_files().read_fits(Pars.full_file, Pars.spall_cols)
df_fits   = read_subset_fits(Pars.dir_fits, Pars.sub_file)
Qsos = Qso_catalog(df_fits)


#test for BOSS and eBOSS
for targ, bits in Qsos.targets.iteritems():
    print "Quasars with only Bit_condition:Ok in {} =".format(targ), Qsos.searching_quasars(targ, bits).sum()


# filter target qsos with Pars.condition
Qsos.filtering_qsos(condition= Pars.condition)


# Compute healpix
unique_pixels = Qsos.adding_pixel_column()


# an example to check whether is working or not
if False: print Qsos.df_qsos.query('THING_ID == 497865723').head()


# If we don't have the files stored, download them from either bnl or sdss website
bnl    = raw_input('Get files from the bnl (y/n): ')
passwd = raw_input('sdss passwd:') if bnl == 'n' else None


# Given a helpix number, find repeated THINGS_ID, and print only those with >= repetitions
print '\n ** {healpix: {THING_ID: num_reps}}'
for _, lpix in enumerate(unique_pixels[:100]):
    thingid_repeat = Qsos.pix_uniqueid(lpix)
    print {lpix: thingid_repeat}
    if not thingid_repeat: continue
    for thids in thingid_repeat:
        #Get specs (from web) given a THING_ID
        qso_files = Qsos.get_files(thing_id= thids, passwd=passwd)
        flag = 1
        while flag:
            print 'doing something'
            #coadd files and compute chisq
            dfall_qsos = Qsos.coadds(qso_files, Pars.spec_cols)
            zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

            #make some plots
            Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
            if flag==1: Qsos.plot_chisq_dist(zipchisq)

            #check specs that have chisq > self.del_chisq, if none, get out
            flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Pars.del_chisq))
            qso_files = Qsos.select_chisq(zipchisq, Pars.del_chisq)
