
"""
Running the qso_catalog.py:

* It goes through the SpAll file and filter all qsos in the BOSS, EBOSS
that satisfy the bit condition and also
    'CLASS== "QSO" & (OBJTYPE=="QSO" | ''OBJTYPE=="NA".) & THING_ID != -1'
* Computes the healpix given the RA and DEC, and group all of them by healpix.
* Withing a healpix, finds the objects with the same THING_ID
* Coadd them (average ivar*flux.) and compute the chisq
* Loop over:  If the chisq is more than 4, eliminate that spec,
    coadd again and get new chisq
"""

from qso_catalog import *
from get_files import *

Pars = Ini_params()

# read the full SpAll file or the subset we're interested on
df_fits = read_subset_fits(Pars.dir_fits, Pars.sub_file)
Qsos    = Qso_catalog(df_fits)

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
Qsos.ask_for_files()

#print_qsos = open('SpAll_files.csv','w')
#print_qsos.write('#Files with repeated THING_ID >= {}\n'.format(Pars.rep_thid))


# Given a helpix number, find repeated THINGS_ID, and print only those with >= repetitions
print '\n ** {healpix: {THING_ID: num_reps}}'
for _, lpix in enumerate(unique_pixels[:100]):
    thingid_repeat = Qsos.pix_uniqueid(lpix)
    if not thingid_repeat: continue
    print {lpix: thingid_repeat}
    for thids in thingid_repeat:
        #Get specs (from web or bnl) given a THING_ID
        qso_files = Qsos.get_files(thing_id= thids)
        #for names in qso_files: print_qsos.write('v5_10_0/spectra/' + names + '\n')
        if False:
            flag = 1
            while flag:
                #coadd files and compute chisq
                dfall_qsos = Qsos.coadds(qso_files, Pars.spec_cols)
                zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                #make some plots
                Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
                if flag==1: Qsos.plot_chisq_dist(zipchisq)

                #check specs that have chisq > self.del_chisq, if none, get out
                flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Pars.del_chisq))
                qso_files = Qsos.select_chisq(zipchisq, Pars.del_chisq)

#print_qsos.close()