
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from get_files import *
from functools import reduce
import fitsio

pd.set_option('display.mpl_style', 'default')

params1 = {'backend': 'pdf',
               'axes.labelsize': 15,
               'text.fontsize': 18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18,
               'legend.fontsize': 8,
               'lines.markersize': 16,
               'font.size': 16,}
#import pylab
#pylab.rcParams.update(params1)


class Ini_params():
    def __init__(self):

        self.trim_chisq  = 2                              #Thresold to discriminate from coadds
        self.rep_thid    = 4                              #Times we want a THING_ID repeated
        self.Npix_side   = 2**5                           #Nside to compute healpix

        self.passwd      = None                           # sdss password
        self.verbose     = False
        self.write_hist  = False                          #Write chisq distribution files
        self.write_names = False                          #Write names of all spec.fits files used
        self.show_plots  = False
        self.use_bokeh   = False                           #Playing with interactive plots

        self.dir_spec    = 'data/spectra/'
        self.dir_v5_10   = 'v5_10_0/spectra/'
        self.pix_dir     = 'healpix/'

        self.full_file   = 'spAll-v5_10_0.fits'
        self.Spall_files = 'SpAll_files.csv'

        self.stats_file  = 'Chisq_dist_all'
        self.stats_file2 = 'Chisq_dist_trim'
        self.stats_file3 = 'Chisq_bad'
        self.suffix      = '_{}.csv'

        self.bit_boss    = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
        self.bit_eboss   = [10,11,12,13,14,15,16,17,18]

        self.targets     = {'BOSS_TARGET1': self.bit_boss, 'EBOSS_TARGET0': self.bit_eboss,
                                   'EBOSS_TARGET1': self.bit_eboss}

        self.condition   = 'CLASS== "QSO".ljust(6) & (OBJTYPE=="QSO".ljust(16) | ' \
                           'OBJTYPE=="NA".ljust(16)) & THING_ID != -1'

        self.spall_cols  = ['RA','DEC','THING_ID','MJD','PLATE','FIBERID','BOSS_TARGET1',
                            'EBOSS_TARGET0','EBOSS_TARGET1']

        self.spec_cols   = ['flux','loglam','ivar','and_mask','or_mask', 'wdisp', 'sky', 'model']

        self.sdss_url    = 'https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/'
        self.bnl_dir     = 'astro:/data/boss/v5_10_0/spectra/'


    def do_nothing(self):
        pass




class Qso_catalog(Ini_params):
    def __init__(self, df_fits, verbose = True):
        self.df_fits     = df_fits
        self.verbose     = verbose
        self.chisq_dist  = []
        Ini_params.__init__(self)



    def searching_quasars(self, data_column, mask_bit):
        """Filter quasars according to the bit array,
        return a Boolean array wherever the quasar is"""

        is_qso  =  lambda bit: ((self.df_fits[data_column] & 2**bit) > 0)
        all_qso =  map(is_qso, mask_bit)
        return reduce(lambda x, y: x | y, all_qso)





    def print_filter_qsos(self, df, text):
        """Print # of Quasars according to filtering condtion"""

        for targ, bits in self.targets.items():
            print ("Quasars with {} condition in {} =".format(text, targ),
                   self.searching_quasars(targ, bits).sum())
        print ("Total qsos = {}".format(len(df)))
        print ("------------------")





    def filtering_qsos(self, condition):
        """Filter quasars within the dataFrame"""

        if self.verbose: self.print_filter_qsos(self.df_fits, 'Bit')

        # only those with CLASS=QSO & OBJTYPE=(QSO|NA)
        self.df_fits  = self.df_fits.query(condition)

        if self.verbose: self.print_filter_qsos(self.df_fits, 'General')

        # and satisfy the bit condition (is a quasar)
        a =[]
        for targ, bits in self.targets.items():
            a.append(self.searching_quasars(targ, bits))
        self.df_qsos = self.df_fits[reduce(lambda x, y: x | y, a)].copy()

        if self.verbose: self.print_filter_qsos(self.df_qsos, 'Both')
        return 0




    def my_own_filter(self, condition):
        """Add your own filter condition"""
        self.df_qsos = self.df_fits.query(condition).copy()




    def adding_pixel_column(self):
        """Computing healpix pixel given 'DEC' and 'RA' """
        phi_rad   = lambda ra : ra*np.pi/180.
        theta_rad = lambda dec: (90.0 - dec)*np.pi/180.

        self.df_qsos['PIX'] = hp.ang2pix(self.Npix_side, theta_rad(self.df_qsos['DEC']), phi_rad(self.df_qsos['RA']))
        unique_pixels  = self.df_qsos['PIX'].unique()
        print ('Unique pixels: ', len(unique_pixels))

        #setting 'PIX' and 'THING_ID' as indices
        self.df_qsos = self.df_qsos.set_index(['PIX','THING_ID'], drop=False).drop('PIX', 1).sort_index()

        if self.verbose: print (self.df_qsos.head())
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
        """nasty hack to get bnl-files, but will change it later"""

        bnl_folder = os.path.join(self.dir_spec, plate)

        if self.verbose: print ('Getting file {} from the bnl'.format(file_name))
        if not os.path.isdir(bnl_folder):
            os.system('mkdir {}'.format(bnl_folder))
        os.system('scp {0} {1}'.format(os.path.join(self.bnl_dir, file_name), bnl_folder))
        return 0





    def get_web_files(self, file_name, passwd):
        """nasty hack to get files from sdss website, but will change it later"""

        if self.verbose: print ('Getting file {} from the sdss web'.format(file_name))
        url = os.path.join(self.sdss_url, file_name)
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

        with open(os.path.join(self.dir_spec, file_name),'wb') as output:
              output.write(data)

        return 0





    def ask_for_files(self, get_files= False):
        """ If we don't have files in local, get them from either bnl
            or sdss website"""

        self.need_files= get_files

        if get_files:
            self.need_files = input('Get files from (bnl), (sdss): ')
            self.passwd     = input('sdss passwd:') if self.need_files == 'sdss' else None

            if not ('{0} == bnl | {0} == sdss'.format(self.need_files)):
                sys.exit('** Need to type either bnl or sdss')
            if self.need_files == 'sdss':
                try:
                    import mechanize
                    from base64 import b64encode
                except:
                    sys.exit("Install mechanize to get files")
        return 0





    def get_names(self, name):
        """ Useful to get values given a THING_ID """
        return list(self.df_qsos.query('THING_ID == %s'%(self.th_id))[name].values)



    def get_files(self, thing_id ='thing_id'):
        """Given a THING_ID locate names of the files,
            if we dont have the files, get them."""

        self.th_id        = thing_id
        self.coadd_id     = 'coadd_%s'%(self.th_id)
        self.ivar_id      = 'ivar_%s'%(self.th_id)
        self.flux_ivar_id = 'flux_ivar_%s'%(self.th_id)

        plates   = self.get_names('PLATE')
        mjds     = self.get_names('MJD')
        fiberids = self.get_names('FIBERID')
        plate_n  = ['{}'.format(plate) for plate in plates]

        qso_files= ['{0}/spec-{0}-{1}-{2}.fits'.format(plate, mjd, str(fiberid).zfill(4))
                        for plate, mjd, fiberid in zip(plates, mjds, fiberids)]

        if self.need_files:
            for plate, file in zip(plate_n, qso_files):
                if not os.path.isfile(self.dir_spec + file):
                    if self.passwd is None:
                        self.get_bnl_files(plate, file)
                    else:
                        self.get_web_files(file, self.passwd)
        return qso_files




    def write_file_names(self):
        """ Write a file that contains all names.fits for filtered quasars,
        so we can get them from NERSC """

        print ('printing names in {}'.format(self.Spall_files))
        self.df_qsos['file_name'] =  self.dir_v5_10 + self.df_qsos['PLATE'].astype(str) + '/spec-' + \
                self.df_qsos['PLATE'].astype(str) + '-' + self.df_qsos['MJD'].astype(str)+ '-' + \
                self.df_qsos['FIBERID'].astype(str).str.zfill(4)

        with open(self.Spall_files, 'w') as f:
             for name in self.df_qsos['file_name'].values:
                 f.write(name + '.fits' + '\n')
        return 0





    def stack_repeated(self, qso_files):
        """Stack all the files in a single one,"""

        stack_qsos = []
        for i, fqso in enumerate(qso_files):
            stack_qsos.append(read_fits(self.dir_spec , fqso, self.spec_cols).set_index('loglam'))
            stack_qsos[i].rename(columns= {'flux': 'flux_%s'%(fqso), 'ivar': 'ivar_%s'%(fqso)}, inplace=True)

        result   = pd.concat([stack_qsos[j][['flux_%s'%(fqso),'ivar_%s'%(fqso)]] for j, fqso in enumerate(qso_files)], axis=1)
        return result.fillna(0).copy()





    def coadds(self, qso_files):
        """ Add coadd column """

        dfall_coadds  = self.stack_repeated(qso_files)
        dfall_coadds[self.flux_ivar_id] = 0
        dfall_coadds[self.ivar_id]      = 0

        for _, fqso in enumerate(qso_files):
            dfall_coadds[self.flux_ivar_id] += dfall_coadds['flux_%s'%(fqso)]*dfall_coadds['ivar_%s'%(fqso)]
            dfall_coadds[self.ivar_id]      += dfall_coadds['ivar_%s'%(fqso)]

        dfall_coadds[self.coadd_id] = dfall_coadds[self.flux_ivar_id]/dfall_coadds[self.ivar_id]

        dfall_coadds = dfall_coadds.fillna(0).copy()
        return dfall_coadds




    def comp_chisq(self, fqso, dfall_coadds):
        """chisq with respect to zero flux """
        tmp = (0.0*dfall_coadds[self.coadd_id].values - dfall_coadds['flux_%s'%(fqso)].values)**2
        return np.sum(tmp*dfall_coadds['ivar_%s'%(fqso)].values)




    def calc_chisq(self, qso_files, dfall_coadds):
        """Compute chisq and return a dict with files'name and chisq"""

        chi_sq_all =[]
        for _, fqso in enumerate(qso_files):
            chis_sq = self.comp_chisq(fqso, dfall_coadds)
            chi_sq_all.append(chis_sq/len(dfall_coadds.values))
                #careful: len(dfall_coadds.values) != len(individuals)
        return dict(zip(qso_files, chi_sq_all))




    def ftrim_chisq(self, zipchisq):
        """Function to select observation that have trim_chisq < chisq"""

        tmp_zipchisq = zipchisq.copy()
        for files, chisq in zipchisq.items():
            if chisq < self.trim_chisq:
                del tmp_zipchisq[files]
            else: continue

        return list(tmp_zipchisq)




    def plot_coadds(self, dfall_coadds, zipchisq):
        """Plot the spectra and coadds"""

        plt.figure(figsize = (18, 8))
        xlimits = [3.55, 4]
        ylimits = [-10, 25]
        ax = plt.subplot(1, 2, 1)
        for fqso, chisq in zipchisq.items():
            dfall_coadds['flux_%s'%(fqso)].plot(label='%s  , Chisq=%s'%(fqso.replace('.fits',''), chisq),
                                         xlim=xlimits, ylim=ylimits, ax=ax)
        plt.legend(loc='best')

        ax2 = plt.subplot(1, 2, 2)
        dfall_coadds[self.coadd_id].plot(label=self.coadd_id, xlim=xlimits, ylim=ylimits, ax=ax2)
        plt.legend(loc='best')
        plt.title('THING_ID: %s'%(self.th_id))
        plt.show(block=True)
        return 0



    def plot_chisq_dist(self, zipchisq):
        """Save chisq for all spectra, and plot a histogram on the fly"""

        for i in zipchisq.values():
            self.chisq_dist.append(i)

        plt.hist(self.chisq_dist, bins=100, range=(0,10))
        plt.ylabel('#')
        plt.xlabel('Chisq')
        plt.title('Chisq Histogram')
        plt.show(block=True)
        #plt.savefig('chisq.pdf')
        return 0



    def write_stats_open(self, rank):
        """Write all chisq and trim after eliminating trim_chisq > chisq"""
        self.write_stats = {'all' : open(self.stats_file  + self.suffix.format(rank), 'w'),
                            'trim': open(self.stats_file2 + self.suffix.format(rank), 'w'),
                            'bad' : open(self.stats_file3 + self.suffix.format(rank), 'w')}


    def write_stats_close(self):
        for i in ['all', 'trim', 'bad']:
            self.write_stats[i].close()


    def write_stats_file(self, zipchisq, name):
        for chi in zipchisq.values():
            self.write_stats[name].write(str(chi) + '\n')



    def write_fits(self, result, lpix):
        data    = (pd.concat([r for r in result], axis=1).fillna(0))
        nrows   = len(data.index)
        data    = data.reset_index().to_dict(orient='list')
        names   = list(data)
        formats = ['f4']*len(names)
        fdata   = np.zeros(nrows, dtype=dict(names= names, formats=formats))
        fdata   = {key:np.array(value) for key,value in data.items()}

        fits = fitsio.FITS(os.path.join(self.pix_dir, 'pix_%s.fits'%(lpix)),'rw')
        fits.write(fdata, header={'Healpix':'%s'%(lpix)})
        fits[1].write_comment("From {}, using Npix_side:{}".format(self.full_file, self.Npix_side))
        if self.verbose: print ('Writing FITS file: %s'%(lpix))
        fits.close()




    def plot_stats(self, size):
        """At the end, collect all chisq and plot a histogram"""
        total_chisq = []
        total_chisq_sec = []
        for i in np.arange(size):
            Chisq     = pd.read_csv(self.stats_file + self.suffix.format(i))
            Chisq_sec = pd.read_csv(self.stats_file2 + self.suffix.format(i))
            total_chisq.extend(Chisq.values.flatten())
            total_chisq_sec.extend(Chisq_sec.values.flatten())

        chisq_ = 'chisq'
        df     = pd.DataFrame(total_chisq,     columns=[chisq_])
        df_sec = pd.DataFrame(total_chisq_sec, columns=[chisq_])

        print (df[chisq_].describe(), df_sec[chisq_].describe())

        bins  = 80
        range = (0,20)

        if self.use_bokeh:
            try:
                from bokeh.plotting import figure, show

                TOOLS = "pan,box_zoom,reset,tap,save,crosshair"
                p1 = figure(title="Counting chisq for repeated THING_ID", tools=TOOLS)

                hist, edges = np.histogram(df[chisq_], density=False, bins=bins, range=range)
                hist2, edges2 = np.histogram(df_sec[chisq_], density=False, bins=bins, range=range)

                p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
                        fill_color="blue", line_color="#FF7373", legend="Chisq-all: %s"%(len(df)))

                p1.quad(top=hist2, bottom=0, left=edges2[:-1], right=edges2[1:],
                        fill_color="red", line_color="#92D7FF", legend="Chisq > 2: %s"%(len(df_sec)))

                p1.xaxis.axis_label = 'Chisq'
                p1.yaxis.axis_label = '#'
                #output_file('histogram.html', title="histogram.py example")
                show(p1)
            except:
                print ('Install Bokeh for fun')
        else:
            plt.figure()
            ax = plt.subplot(111)
            df[chisq_].plot.hist(    bins=bins, range=range, alpha=0.9, ax=ax, color='r', label='Chisqs, %s'%(len(df)))
            df_sec[chisq_].plot.hist(bins=bins, range=range, alpha=0.5, ax=ax, color='b', label='After loop, %s'%(len(df_sec)))

            plt.ylabel('#')
            plt.xlabel(chisq_)
            plt.legend(loc = 'best')

            plt.show(block=True)




if __name__=='__main__':
    print ("goofing around :P ")



