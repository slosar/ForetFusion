
import os, sys
import pylab
import numpy as np
import healpy as hp
import pandas as pd
import Get_files
import matplotlib.pyplot as plt

pd.set_option('display.mpl_style', 'default')

params1 = {'backend': 'pdf',
               'axes.labelsize': 15,
               'text.fontsize': 18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18,
               'legend.fontsize': 8,
               'lines.markersize': 16,
               'font.size': 16,}
pylab.rcParams.update(params1)


class Qso_catalog():
    def __init__(self, df_fits, verbose = True):
        self.df_fits     = df_fits
        self.verbose     = verbose
        self.Npix_side   = 2**5


    def searching_quasars(self, data_column, mask_bit):
        """Filter the quasar according to the bit array"""

        is_qso  =  lambda bit: ((self.df_fits[data_column] & 2**bit) > 0)
        all_qso =  map(is_qso, mask_bit)
        return reduce(lambda x, y: x | y, all_qso)



    def filtering_qsos(self, targets, condition):
        """Filter only the quasars"""

        # only those with CLASS=QSO & OBJTYPE=(QSO|NA)
        self.df_fits  = self.df_fits.query(condition)

        # and satisfy the bit condition
        a =[]
        for targ, bits in targets.iteritems():
            a.append(self.searching_quasars(targ, bits))
        self.df_qsos = self.df_fits[reduce(lambda x, y: x | y, a)].copy()


        print 'Total qsos (Bit_condition:Ok & %s='%(condition), len(self.df_qsos)
        return 0




    def adding_pixel_column(self):
        """Computing healpix pixel given 'DEC' and 'RA' """
        phi_rad   = lambda ra : ra*np.pi/180.
        theta_rad = lambda dec: (90.0 - dec)*np.pi/180.

        self.df_qsos['PIX'] = hp.ang2pix(self.Npix_side, theta_rad(self.df_qsos['DEC']), phi_rad(self.df_qsos['RA']))
        self.unique_pixels  = self.df_qsos['PIX'].unique()
        print 'Unique pixels: ', len(self.unique_pixels)

        #setting 'PIX' and 'THING_ID' as indices
        self.df_qsos = self.df_qsos.set_index(['PIX','THING_ID'], drop=False).drop('PIX', 1).sort_index()

        if self.verbose: print self.df_qsos.head()
        return 0




    def pix_uniqueid(self, pix_id, repetitions= 2):
        """Given a pixel, return the THING_ID and the number of times is repeated"""
        rep_thing_id = self.df_qsos.query('PIX == {}'.format(pix_id)).groupby('THING_ID').size()

        if not rep_thing_id[rep_thing_id >= repetitions].empty:
            self.uniqeid = dict(rep_thing_id[rep_thing_id >= repetitions])
        else:
            self.uniqeid = 0
        return self.uniqeid




    def get_names(self, thing_id, name):
        return list(self.df_qsos.query('THING_ID == %s'%(thing_id))[name].values)




    def get_files(self, thing_id ='thing_id'):
        plates   = self.get_names(thing_id, 'PLATE')
        mjds     = self.get_names(thing_id, 'MJD')
        fiberids = self.get_names(thing_id, 'FIBERID')
        plate_n  = ['%s'%(plate) for plate in plates]

        files    = ['spec-%s-%s-%s.fits'%(plate, mjd, str(fiberid).zfill(4))
                        for plate, mjd, fiberid in zip(plates, mjds, fiberids)]


        web_file = []
        for plate, file in zip(plate_n, files):
            print 'Getting file {} from the web'.format(file)
            if not os.path.isfile(file):
                Get_files.get_web_files(plate, file)

            web_file.append('v5_10_0/spectra/%s/%s'%(plate, file))

        return files, web_file



    def stack_repeated(self, qsos_files, columns):
        stack_qsos = []
        for i, stacks in enumerate(qsos_files):
            stack_qsos.append(Get_files.read_fits(stacks, columns).set_index('loglam'))
            stack_qsos[i]['flux_%s'%(i)] = stack_qsos[i]['flux']
            stack_qsos[i]['ivar_%s'%(i)] = stack_qsos[i]['ivar']

        result   = pd.concat([stack_qsos[j][['flux_%s'%(j),'ivar_%s'%(j)]] for j, _ in enumerate(qsos_files)], axis=1)
        return result.fillna(0).copy()



    def coadds(self, qsos_files, columns):
        all_qsos = self.stack_repeated(qsos_files, columns)
        all_qsos['sum_flux_ivar']=0
        all_qsos['sum_ivar']=0
        for i, _ in enumerate(qsos_files):
            all_qsos['sum_flux_ivar'] += all_qsos['flux_%s'%(i)]*all_qsos['ivar_%s'%(i)]
            all_qsos['sum_ivar']      += all_qsos['ivar_%s'%(i)]
        all_qsos['coadd'] = all_qsos['sum_flux_ivar']/all_qsos['sum_ivar']

        return all_qsos


    def plot_coadds(self, thingid, all_qsos):
        plt.figure(figsize = (18, 8))
        xlimits = [3.55,4]
        ylimits = [-10,20]
        ax = plt.subplot(1,2,1)
        for i,_ in enumerate(qsos_files):
            all_qsos['flux_%s'%(i)].plot(label=qsos_files[i],xlim=xlimits, ylim=ylimits, ax=ax)
        plt.legend(loc='best')

        ax2 = plt.subplot(1,2,2)
        all_qsos['coadd'].plot(label='coad', xlim=xlimits, ylim=ylimits, ax=ax2)
        plt.legend(loc='best')
        plt.title('THING_ID: %s'%(thingid))
        plt.show(block=True)


if __name__=='__main__':

    #read the full file
    if False:
        file_name = 'spAll-v5_10_0.fits'
        columns   = ['RA','DEC','THING_ID','MJD','PLATE','FIBERID','BOSS_TARGET1','EBOSS_TARGET0','EBOSS_TARGET1', 'etc....']
        df_qsos   = Get_files.read_fits(file_name, columns)

    #instead read the subset we're interested on
    df_fits   = Get_files.read_subset_fits('subset_spAll-v5_10_0.csv')

    Qsos      = Qso_catalog(df_fits, verbose= False)

    bit_boss  = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
    bit_eboss = [10,11,12,13,14,15,16,17,18]
    targets   = {'BOSS_TARGET1': bit_boss, 'EBOSS_TARGET0': bit_eboss, 'EBOSS_TARGET1': bit_eboss}
    for targ, bits in targets.iteritems():
        print 'Quasars with Bit_condition:Ok in {} ='.format(targ), Qsos.searching_quasars(targ, bits).sum()


    # filter qsos with CLASS,OBJTYPE and THING_ID != -1
    condition = 'CLASS== "QSO".ljust(6) & (OBJTYPE=="QSO".ljust(16) | OBJTYPE=="NA".ljust(16)) & THING_ID != -1'
    Qsos.filtering_qsos(targets, condition= condition)

    # Compute healpix
    Qsos.adding_pixel_column()


    if False:
        print '\n  *** just checking'
        print Qsos.df_qsos.query('THING_ID == 497865723').head()

    spec_columns    = ['flux','loglam','ivar','and_mask','or_mask', 'sky', 'wdisp', 'model']

    #f = open('SpAll_files.csv','w')
    #f.write('#Files with repeated THING_ID \n')
    # Given a pixel number, find repeated THINGS_ID, and print only >= reps
    for i, lpix in enumerate(Qsos.unique_pixels):
        thingid_repeat = Qsos.pix_uniqueid(lpix, repetitions= 2)
        print {lpix: thingid_repeat}
        #print i, lpix
        if thingid_repeat != 0:
            for ids in thingid_repeat:
                 ## Get files (from web) given a THING_ID
                qsos_files, web_file = Qsos.get_files(thing_id= ids)
                all_qsos   = Qsos.coadds(qsos_files, spec_columns)
                Qsos.plot_coadds(ids, all_qsos)

                #for i in web_file:
                 #   f.write(i+'\n')
