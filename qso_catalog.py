
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
        self.del_chisq   = 4

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

        qso_files= ['spec-%s-%s-%s'%(plate, mjd, str(fiberid).zfill(4))
                        for plate, mjd, fiberid in zip(plates, mjds, fiberids)]

        for plate, file in zip(plate_n, qso_files):
            file += '.fits'
            if not os.path.isfile(file):
                print 'Getting file {} from the web'.format(file)
                Get_files.get_web_files(plate, file)

        return qso_files



    def stack_repeated(self, qso_files, columns):
        stack_qsos = []
        for i, fqso in enumerate(qso_files):
            stack_qsos.append(Get_files.read_fits(fqso, columns).set_index('loglam'))
            stack_qsos[i]['flux_%s'%(fqso)] = stack_qsos[i]['flux']
            stack_qsos[i]['ivar_%s'%(fqso)] = stack_qsos[i]['ivar']

        result   = pd.concat([stack_qsos[j][['flux_%s'%(stacks),'ivar_%s'%(stacks)]] for j, stacks in enumerate(qso_files)], axis=1)
        return result.fillna(0).copy()



    def coadds(self, qso_files, columns):
        dfall_qsos = self.stack_repeated(qso_files, columns)
        dfall_qsos['sum_flux_ivar'] =0
        dfall_qsos['sum_ivar']      =0
        for i, fqso in enumerate(qso_files):
            dfall_qsos['sum_flux_ivar'] += dfall_qsos['flux_%s'%(fqso)]*dfall_qsos['ivar_%s'%(fqso)]
            dfall_qsos['sum_ivar']      += dfall_qsos['ivar_%s'%(fqso)]
        dfall_qsos['coadd'] = dfall_qsos['sum_flux_ivar']/dfall_qsos['sum_ivar']

        dfall_qsos = dfall_qsos.fillna(0).copy()
        return dfall_qsos



    def calc_chisq(self, qso_files, dfall_qsos):
        chi_sq_all =[]
        for i, fqso in enumerate(qso_files):
            chis_sq = np.sum((dfall_qsos['coadd'].values - dfall_qsos['flux_%s'%(fqso)].values)**2*dfall_qsos['ivar_%s'%(fqso)].values)
            chi_sq_all.append(chis_sq/len(dfall_qsos.values))
        return dict(zip(qso_files, chi_sq_all))



    def select_chisq(self, zipchisq):
        rm_file = []
        for files, chisq in zipchisq.iteritems():
            if chisq > self.del_chisq:
                rm_file.append(files)

        for rm in rm_file:
            del zipchisq[rm]

        return zipchisq.keys()



    def plot_coadds(self, dfall_qsos, thingid, zipchisq):
        plt.figure(figsize = (18, 8))
        xlimits = [3.55, 4]
        ylimits = [-10, 25]
        ax = plt.subplot(1, 2, 1)
        for fqso, chisq in zipchisq.iteritems():
            dfall_qsos['flux_%s'%(fqso)].plot(label='%s  , chisq=%s'%(fqso, chisq),
                                         xlim=xlimits, ylim=ylimits, ax=ax)
        plt.legend(loc='best')

        ax2 = plt.subplot(1,2,2)
        dfall_qsos['coadd'].plot(label='coad', xlim=xlimits, ylim=ylimits, ax=ax2)
        plt.legend(loc='best')
        plt.title('THING_ID: %s'%(thingid))
        plt.show(block=True)
        return 0




class Ini_params():
    def __init__(self):

        self.verbose   = False
        self.sub_file  = 'subset_spAll-v5_10_0.csv'
        self.full_file = 'spAll-v5_10_0.fits'

        self.bit_boss  = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
        self.bit_eboss = [10,11,12,13,14,15,16,17,18]

        self.condition = 'CLASS== "QSO".ljust(6) & (OBJTYPE=="QSO".ljust(16) | ' \
                         'OBJTYPE=="NA".ljust(16)) & THING_ID != -1'

        self.targets   = {'BOSS_TARGET1': self.bit_boss, 'EBOSS_TARGET0': self.bit_eboss,
                                   'EBOSS_TARGET1': self.bit_eboss}

        self.spall_cols = ['RA','DEC','THING_ID','MJD','PLATE','FIBERID','BOSS_TARGET1',
                           'EBOSS_TARGET0','EBOSS_TARGET1']

        self.spec_cols  = ['flux','loglam','ivar','and_mask','or_mask', 'sky', 'wdisp', 'model']


    def do_nothing(self):
        pass


if __name__=='__main__':
    Pars      = Ini_params()

    #read the full file
    if False: df_qsos   = Get_files.read_fits(Pars.full_file, Pars.spall_cols)

    #instead read the subset we're interested on
    df_fits   = Get_files.read_subset_fits(Pars.sub_file)
    Qsos      = Qso_catalog(df_fits, verbose= Pars.verbose)

    for targ, bits in Pars.targets.iteritems():
        print 'Quasars with Bit_condition:Ok in {} ='.format(targ), Qsos.searching_quasars(targ, bits).sum()

    # filter qsos with CLASS,OBJTYPE and THING_ID != -1
    Qsos.filtering_qsos(Pars.targets, condition= Pars.condition)

    # Compute healpix
    Qsos.adding_pixel_column()

    if False: print Qsos.df_qsos.query('THING_ID == 497865723').head()

    # Given a helpix number, find repeated THINGS_ID, and print only those with >= repetitions
    print '\n ** {healpix: {THING_ID: num_reps}}'
    for _, lpix in enumerate(Qsos.unique_pixels):
        thingid_repeat = Qsos.pix_uniqueid(lpix, repetitions= 5)
        print {lpix: thingid_repeat}
        if thingid_repeat != 0:
            for thids in thingid_repeat:
                flag = 1
                #Get specs (from web) given a THING_ID
                qso_files = Qsos.get_files(thing_id= thids)
                while flag:
                    #coadd files and compute chisq
                    dfall_qsos = Qsos.coadds(qso_files, Pars.spec_cols)
                    zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                    #make some plots
                    Qsos.plot_coadds(dfall_qsos, thids, zipchisq)

                    #check specs that have chisq > self.del_chisq
                    #if none, get out
                    flag = len(qso_files) - len(Qsos.select_chisq(zipchisq))
                    qso_files = Qsos.select_chisq(zipchisq)


