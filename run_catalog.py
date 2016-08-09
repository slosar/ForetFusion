
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
import Get_files

# import initial parameters
Pars = Ini_params()

# read the full SpAll file
if False: df_qsos = Get_files.read_fits(Pars.direct, Pars.full_file, Pars.spall_cols)

# instead read the subset we're interested on
df_fits   = Get_files.read_subset_fits(Pars.direct, Pars.sub_file)
Qsos      = Qso_catalog(df_fits, verbose= Pars.verbose)


for targ, bits in Pars.targets.iteritems():
    print "Quasars with Bit_condition:Ok in {} =".format(targ), Qsos.searching_quasars(targ, bits).sum()

# filter target qsos with Pars.condition
Qsos.filtering_qsos(Pars.targets, condition= Pars.condition)

# Compute healpix
unique_pixels = Qsos.adding_pixel_column(Pars.Npix_side)

# an example to check whether is working or not
if False: print Qsos.df_qsos.query('THING_ID == 497865723').head()

# Given a helpix number, find repeated THINGS_ID, and print only those with >= repetitions
print '\n ** {healpix: {THING_ID: num_reps}}'

bnl    = raw_input('Are you @ bnl (y/n): ')
passwd = raw_input('sdss passwd:') if bnl == 'n' else None


for _, lpix in enumerate(unique_pixels[:100]):
    thingid_repeat = Qsos.pix_uniqueid(lpix, repetitions= Pars.rep_thid)
    print {lpix: thingid_repeat}
    if not thingid_repeat: continue
    for thids in thingid_repeat:
        #Get specs (from web) given a THING_ID
        qso_files, plate = Qsos.get_files(Pars.direct, thing_id= thids, passwd=passwd)
        flag = 1
        while flag:
            print 'doing something'
            #coadd files and compute chisq
            dfall_qsos = Qsos.coadds(Pars.direct, plate, qso_files, Pars.spec_cols)
            zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

            #make some plots
            Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
            if flag==1: Qsos.plot_chisq_dist(zipchisq)

            #check specs that have chisq > self.del_chisq, if none, get out
            flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Pars.del_chisq))
            qso_files = Qsos.select_chisq(zipchisq, Pars.del_chisq)
