
import pylab
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import mechanize
from base64 import b64encode
from get_files import *

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


class Ini_params():
    def __init__(self):

        self.del_chisq  = 4                              #Thresold to discriminate from coadds
        self.rep_thid   = 4                              #Times we want a THING_ID repeated
        self.Npix_side  = 2**5                           #Nside to compute healpix
        self.need_files = 'No'                            #Assume we have the files
        self.passwd     = None                           # sdss password
        self.verbose    = False
        self.print_stat = True
        self.dir_fits   = 'data/'
        self.dir_spec   = self.dir_fits + 'spectra/'
        self.full_file  = 'spAll-v5_10_0.fits'
        self.sub_file   = 'subset_spAll-v5_10_0.csv'

        self.bit_boss  = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
        self.bit_eboss = [10,11,12,13,14,15,16,17,18]

        self.targets   = {'BOSS_TARGET1': self.bit_boss, 'EBOSS_TARGET0': self.bit_eboss,
                                   'EBOSS_TARGET1': self.bit_eboss}

        self.condition = 'CLASS== "QSO".ljust(6) & (OBJTYPE=="QSO".ljust(16) | ' \
                         'OBJTYPE=="NA".ljust(16)) & THING_ID != -1'

        self.spall_cols = ['RA','DEC','THING_ID','MJD','PLATE','FIBERID','BOSS_TARGET1',
                           'EBOSS_TARGET0','EBOSS_TARGET1']

        self.spec_cols  = ['flux','loglam','ivar','and_mask','or_mask', 'wdisp', 'sky', 'model']

    def do_nothing(self):
        pass





class Qso_catalog(Ini_params):
    def __init__(self, df_fits, verbose = True):
        self.df_fits     = df_fits
        self.verbose     = verbose
        self.chisq_dist  = []
        Ini_params.__init__(self)



    def searching_quasars(self, data_column, mask_bit):
        """Filter the quasar according to the bit array,
        return a Boolean array wherever the quasar is"""

        is_qso  =  lambda bit: ((self.df_fits[data_column] & 2**bit) > 0)
        all_qso =  map(is_qso, mask_bit)
        return reduce(lambda x, y: x | y, all_qso)



    def filtering_qsos(self, condition):
        """Filter only the quasars withing the dataFrame"""

        # only those with CLASS=QSO & OBJTYPE=(QSO|NA)
        self.df_fits  = self.df_fits.query(condition)

        # and satisfy the bit condition
        a =[]
        for targ, bits in self.targets.iteritems():
            a.append(self.searching_quasars(targ, bits))
        self.df_qsos = self.df_fits[reduce(lambda x, y: x | y, a)].copy()


        print 'Total qsos (Bit_condition:Ok & %s='%(condition), len(self.df_qsos)
        return 0




    def adding_pixel_column(self):
        """Computing healpix pixel given 'DEC' and 'RA' """
        phi_rad   = lambda ra : ra*np.pi/180.
        theta_rad = lambda dec: (90.0 - dec)*np.pi/180.

        self.df_qsos['PIX'] = hp.ang2pix(self.Npix_side, theta_rad(self.df_qsos['DEC']), phi_rad(self.df_qsos['RA']))
        unique_pixels  = self.df_qsos['PIX'].unique()
        print 'Unique pixels: ', len(unique_pixels)

        #setting 'PIX' and 'THING_ID' as indices
        self.df_qsos = self.df_qsos.set_index(['PIX','THING_ID'], drop=False).drop('PIX', 1).sort_index()

        if self.verbose: print self.df_qsos.head()
        return unique_pixels



    def pix_uniqueid(self, pix_id):
        """Given a pixel, return the THING_ID and the number of times is repeated"""
        rep_thing_id = self.df_qsos.query('PIX == {}'.format(pix_id)).groupby('THING_ID').size()

        if not rep_thing_id[rep_thing_id >= self.rep_thid].empty:
            uniqeid = dict(rep_thing_id[rep_thing_id >= self.rep_thid])
        else:
            uniqeid = {}
        return uniqeid



    def get_bnl_files(self, plate, file_name):
        """nasty hack, but change it later"""
        print 'Getting file {} from the bnl'.format(file_name)
        if not os.path.isdir('{}{}'.format(self.dir_spec, plate)):
            os.system('mkdir {}{}'.format(self.dir_spec, plate))
        os.system('scp astro:/data/boss/v5_10_0/spectra/{1} {2}{0}'.format(plate, file_name, self.dir_spec))
        return 0



    def get_web_files(self, file_name, passwd):
        print 'Getting file {} from the web'.format(file_name)
        url = 'https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/{}'.format(file_name)
        username = 'sdss'
        password = '{}'.format(passwd)

        # I have had to add a carriage return ('%s:%s\n'), but
        # you may not have to.
        b64login = b64encode('%s:%s' % (username, password))

        br = mechanize.Browser()
        br.set_handle_robots(False)

        br.addheaders.append(
          ('Authorization', 'Basic %s' % b64login )
        )
        br.open(url)
        r = br.response()
        data = r.read()

        with open('{}{}'.format(self.dir_spec, file_name),'wb') as output:
              output.write(data)

        return 0



    def ask_for_files(self, get_them=False):
        self.need_files= get_them
        if get_them:
            self.need_files = raw_input('Get files from bnl(1), sdss(2): ')
            self.passwd     = raw_input('sdss passwd:') if self.need_files == '2' else None

        return 0


    def get_names(self, thing_id, name):
        return list(self.df_qsos.query('THING_ID == %s'%(thing_id))[name].values)



    def get_files(self, thing_id ='thing_id'):
        plates   = self.get_names(thing_id, 'PLATE')
        mjds     = self.get_names(thing_id, 'MJD')
        fiberids = self.get_names(thing_id, 'FIBERID')
        plate_n  = ['{}'.format(plate) for plate in plates]

        qso_files= ['{0}/spec-{0}-{1}-{2}'.format(plate, mjd, str(fiberid).zfill(4))
                        for plate, mjd, fiberid in zip(plates, mjds, fiberids)]

        if self.need_files:
            for plate, file in zip(plate_n, qso_files):
                file = '{}.fits'.format(file)
                if not os.path.isfile(self.dir_spec + file):
                    if self.passwd is None:
                        self.get_bnl_files(plate, file)
                    else:
                        self.get_web_files(file, self.passwd)
        return qso_files


    def print_file_names(self):
        nfile = 'SpAll_files.csv'
        print 'printing names in {}'.format(nfile)
        self.df_qsos['file_name'] =  'v5_10_0/spectra/' + self.df_qsos['PLATE'].astype(str) + '/spec-' + \
                                     self.df_qsos['PLATE'].astype(str) + '-' + self.df_qsos['MJD'].astype(str)+ '-' + \
                                     self.df_qsos['FIBERID'].astype(str).str.zfill(4)

        with open(nfile ,'w') as f:
             for name in self.df_qsos['file_name'].values:
                 f.write(name + '\n')
        return 0



    def stack_repeated(self, qso_files, columns):
        stack_qsos = []
        for i, fqso in enumerate(qso_files):
            stack_qsos.append(read_fits(self.dir_spec, fqso, columns).set_index('loglam'))
            stack_qsos[i]['flux_%s'%(fqso)] = stack_qsos[i]['flux']
            stack_qsos[i]['ivar_%s'%(fqso)] = stack_qsos[i]['ivar']

        result   = pd.concat([stack_qsos[j][['flux_%s'%(stacks),'ivar_%s'%(stacks)]] for j, stacks in enumerate(qso_files)], axis=1)
        return result.fillna(0).copy()



    def coadds(self, qso_files, columns):
        dfall_qsos = self.stack_repeated(qso_files, columns)
        dfall_qsos['sum_flux_ivar'] = 0
        dfall_qsos['sum_ivar']      = 0
        for i, fqso in enumerate(qso_files):
            dfall_qsos['sum_flux_ivar'] += dfall_qsos['flux_%s'%(fqso)]*dfall_qsos['ivar_%s'%(fqso)]
            dfall_qsos['sum_ivar']      += dfall_qsos['ivar_%s'%(fqso)]

        dfall_qsos['coadd'] = dfall_qsos['sum_flux_ivar']/dfall_qsos['sum_ivar']

        dfall_qsos = dfall_qsos.fillna(0).copy()
        return dfall_qsos



    def calc_chisq(self, qso_files, dfall_qsos):
        """Compute chisq and return a dict with files'name and chisq"""
        chi_sq_all =[]
        for i, fqso in enumerate(qso_files):
            chis_sq = np.sum((dfall_qsos['coadd'].values - dfall_qsos['flux_%s'%(fqso)].values)**2*dfall_qsos['ivar_%s'%(fqso)].values)
            chi_sq_all.append(chis_sq/len(dfall_qsos.values))
        return dict(zip(qso_files, chi_sq_all))



    def select_chisq(self, zipchisq, del_chisq):
        tmp_zipchisq = zipchisq.copy()
        for files, chisq in zipchisq.iteritems():
            if chisq > del_chisq:
                del tmp_zipchisq[files]
            else: continue

        return tmp_zipchisq.keys()



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


    def plot_chisq_dist(self, zipchisq):
        for i in zipchisq.values():
            self.chisq_dist.append(i)

        plt.hist(self.chisq_dist, bins=100, range=(0,10))
        plt.ylabel('#')
        plt.xlabel('chisq')
        plt.title('chisq Histogram')
        #plt.show(block=True)
        plt.savefig('chisq.pdf')
        return 0



if __name__=='__main__':
    Pars      = Ini_params()


    #instead read the subset we're interested on
    df_fits   = Get_files.read_subset_fits(Pars.sub_file)
    Qsos      = Qso_catalog(df_fits, verbose= Pars.verbose)

    for targ, bits in Pars.targets.iteritems():
        print 'Quasars with Bit_condition:Ok in {} ='.format(targ), Qsos.searching_quasars(targ, bits).sum()



