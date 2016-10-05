
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from get_files import *
from functools import reduce
import seaborn as sns
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
        self.Npix_side   = 2**6                           #Nside to compute healpix

        self.passwd      = None                           #sdss password
        self.run_sky     = False                          #to run sky+flux calculations
        self.write_master= False                          #Write master file with info
        self.write_ffits = False                          #Write fits files for each pix
        self.write_hist  = False                          #Write chisq distribution files
        self.write_names = False                          #Write names of all spec.fits files used
        self.show_plots  = False                          #must be False when using mpi
        self.use_bokeh   = False                          #Playing with interactive plots

        self.dir_spec    = 'data/spectra/'
        self.dir_v5_10   = 'v5_10_0/spectra/'
        self.pix_dir     = 'healpix/'

        self.full_file   = 'spAll-v5_10_0.fits'
        self.Spall_files = 'SpAll_files.csv'

        self.stats_file  = 'Chisq_dist'
        self.stats_file3 = 'Chisq_bad'
        self.suffix      = '_{}.csv'

        self.bit_boss    = [10,11,12,13,14,15,16,17,18,19,40,41,42,43,44]
        self.bit_eboss   = [10,11,12,13,14,15,16,17,18]

        self.targets     = {'BOSS_TARGET1': self.bit_boss, 'EBOSS_TARGET0': self.bit_eboss,
                                   'EBOSS_TARGET1': self.bit_eboss}

        self.condition   = 'CLASS== "QSO".ljust(6) & (OBJTYPE=="QSO".ljust(16) | ' \
                           'OBJTYPE=="NA".ljust(16)) & THING_ID != -1'

        self.spall_cols  = ['RA','DEC','THING_ID','MJD','PLATE','FIBERID','BOSS_TARGET1','CLASS','OBJTYPE',
                            'EBOSS_TARGET0','EBOSS_TARGET1','Z','Z_ERR','ZWARNING','MODELMAG']

        self.spec_cols   = ['flux','loglam','ivar','and_mask','or_mask', 'wdisp', 'sky', 'model']

        self.sdss_url    = 'https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/'
        self.bnl_dir     = 'astro:/data/boss/v5_10_0/spectra/'


    def do_nothing(self):
        pass




class Qso_catalog(Ini_params):
    def __init__(self, df_fits, verbose = True):
        self.df_fits     = df_fits
        self.verbose     = verbose
        self.df_qsos     = None         #Main DataFrame that contains all Qso info.
        self.chisq_dist  = []
        self.all_lpix    = []
        self.all_thid    = []
        self.all_qfiles  = []
        Ini_params.__init__(self)



    def searching_quasars(self, data_column, mask_bit):
        """Filter quasars according to the bit array,
        return a Boolean array wherever the quasar is"""

        is_qso  =  lambda bit: ((self.df_fits[data_column] & 2**bit) > 0)
        all_qso =  map(is_qso, mask_bit)
        return reduce(lambda x, y: x | y, all_qso)





    def print_filter_qsos(self, df, text):
        """Print # of Quasars according to the filtering condtion"""

        for targ, bits in self.targets.items():
            print ("Quasars with {} condition in {} =".format(text, targ),
                   self.searching_quasars(targ, bits).sum())
        print ("Total qsos = {}".format(len(df)))
        print ("------------------")





    def filtering_qsos(self, condition):
        """Filter quasars given the condition"""

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
        return 0




    def adding_pixel_column(self):
        """Computing healpix pixel given 'DEC' and 'RA', and return unique healpix """

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
        """nasty hack to get files from sdss website, but will change it later,
            only works with Python2"""

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
        """ If we don't have files in local directory, get them from either bnl
            or sdss website"""

        self.need_files= get_files

        if get_files:
            self.need_files = input('Get files from (bnl), (sdss): ')
            self.passwd     = input('sdss passwd:') if self.need_files == 'sdss' else None

            if not ('{0} == bnl | {0} == sdss'.format(self.need_files)):
                os.sys.exit('** Need to type either bnl or sdss')
            if self.need_files == 'sdss':
                try:
                    import mechanize
                    from base64 import b64encode
                except:
                    os.sys.exit("Install mechanize to get files")
        return 0




    def get_names(self, name):
        """ Useful to get values given a THING_ID """
        return list(self.df_qsos.query('THING_ID == %s'%(self.th_id))[name].values)




    def get_files(self, thing_id ='thing_id'):
        """Given a THING_ID locate names of the files,
            if we dont have the files, get them."""

        self.th_id        = thing_id

        plates   = self.get_names('PLATE')
        mjds     = self.get_names('MJD')
        fiberids = self.get_names('FIBERID')
        plate_nu = ['{}'.format(plate) for plate in plates]

        qso_files= ['{0}/spec-{0}-{1}-{2}.fits'.format(plate, mjd, str(fiberid).zfill(4))
                        for plate, mjd, fiberid in zip(plates, mjds, fiberids)]

        #just in case we need to have the files stored
        if self.need_files:
            for plate, file in zip(plate_nu, qso_files):
                if not os.path.isfile(self.dir_spec + file):
                    if self.passwd is None:
                        self.get_bnl_files(plate, file)
                    else:
                        self.get_web_files(file, self.passwd)
        return qso_files




#    def get_zvalues(self, thing_id ='thing_id'):
#	""" Get ZWARNING, Z_ERR, Z"""
	
#	zwarning = self.get_names('ZWARNING')
#	zerr	 = self.get_names('Z_ERR')	
#	z	 = self.get_names('Z')
#	return (zwarning, zerr, z)




    def write_file_names(self):
        """ Write a file that contains all names.fits for filtered quasars,
        so we can get them later from NERSC """

        fnames = 'fnames'

        print ('... Printing file-names in {}'.format(self.Spall_files))

        self.df_qsos[fnames] =  self.df_qsos['PLATE'].astype(str) + '/spec-' + \
                self.df_qsos['PLATE'].astype(str) + '-' + self.df_qsos['MJD'].astype(str)+ '-' + \
                self.df_qsos['FIBERID'].astype(str).str.zfill(4) + '.fits'

        with open(self.Spall_files, 'w') as f:
             for name in (self.dir_v5_10 + self.df_qsos[fnames]).values:
                 f.write(name + '\n')

        return self.df_qsos[fnames]





    def stack_them(self, dic_file):
        """Stack all the files in a single one, so is easier to coadd them"""

        stack_qsos = []
        for fqso, df in dic_file.items():
            columns= {'flux': 'flux_%s'%(fqso), 'ivar': 'ivar_%s'%(fqso), 'sky': 'sky_%s'%(fqso),
                        'and_mask': 'and_mask_%s'%(fqso), 'or_mask': 'or_mask_%s'%(fqso)}
            df.rename(columns= columns, inplace=True)
            stack_qsos.append(df)

        result = pd.concat(stack_qsos, axis=1)
        return result.fillna(0).copy()





    def coadds(self, dic_file):
        """ Add coadd, and_mask and or_mask columns for a given THING_ID """

        self.coadd_id     = 'coadd_%s'%(self.th_id)
        self.ivar_id      = 'ivar_%s'%(self.th_id)
        self.flux_ivar_id = 'flux_ivar_%s'%(self.th_id)
        self.and_mask_id  = 'and_mask_%s'%(self.th_id)
        self.or_mask_id   = 'or_mask_%s'%(self.th_id)


        dfall_coadds  = self.stack_them(dic_file)
        dfall_coadds[self.flux_ivar_id] = 0
        dfall_coadds[self.ivar_id]      = 0

        for _, fqso in enumerate(dic_file.keys()):
            flux = dfall_coadds['flux_%s'%(fqso)]
            if self.run_sky: flux += dfall_coadds['sky_%s'%(fqso)]
            dfall_coadds[self.flux_ivar_id]    += ( flux  )*dfall_coadds['ivar_%s'%(fqso)]
            dfall_coadds[self.ivar_id]         += dfall_coadds['ivar_%s'%(fqso)]
            dfall_coadds['or_mask_%s'%(fqso)]   = pd.DataFrame(dfall_coadds['or_mask_%s'%(fqso)], dtype='int')
            dfall_coadds['and_mask_%s'%(fqso)]  = pd.DataFrame(dfall_coadds['and_mask_%s'%(fqso)], dtype='int')


        dfall_coadds[self.and_mask_id] = (reduce(lambda x, y: x & y, [dfall_coadds['and_mask_%s'%(i)] for i in dic_file.keys()]))
        dfall_coadds[self.or_mask_id]  = (reduce(lambda x, y: x | y, [dfall_coadds['or_mask_%s'%(i)]  for i in dic_file.keys()]))
        dfall_coadds[self.coadd_id]    = dfall_coadds[self.flux_ivar_id] / dfall_coadds[self.ivar_id]

        dfall_coadds = dfall_coadds.fillna(0).copy()
        return dfall_coadds



    def comp_chisq(self, df):
        return np.sum((df['flux'].values)**2*df['ivar'].values)/len(df['flux'].values)



    def cal_chisq(self, qso_files):
        dic_file  = {}
        dic_chisq = {}
	dict_z    = self.get_zvalues(qso_files)

        for fqso in qso_files:
            df_file = read_fits(self.dir_spec , fqso, self.spec_cols).set_index('loglam')
            chisq   = self.comp_chisq(df_file)
            if chisq > self.trim_chisq:
                dic_file[fqso]  = df_file
                dic_chisq[fqso] = chisq
	    else:
		del dict_z[fqso]
	        if self.write_hist: self.write_stats_file(str(self.th_id) + ',' + str(fqso) + ', chisq: ' + str(chisq), 'bad')
        return dic_file, dic_chisq, dict_z



    def get_zvalues(self, qso_files):
        """ Get ZWARNING, Z_ERR, Z"""
	dict_z = {}
        zwarning = self.get_names('ZWARNING')
        zerr     = self.get_names('Z_ERR')
        z        = self.get_names('Z')
	for file, zwar, zerr, z in zip(qso_files, zwarning, zerr, z):
	    dict_z[file]=(zwar, zerr, z)	
        return dict_z





    def plot_coadds(self, dfall_qsos, dic_chisq):
        """Plot the spectra and coadds"""

        plt.figure(figsize = (18, 8))
        xlimits = [3.55, 4]; ylimits = [-10, 25]

        ax = plt.subplot(1, 2, 1)
        for fqso in dic_chisq.keys():
            flux = 'flux_%s'%(fqso)
            if self.run_sky: flux = [flux, 'sky_%s'%(fqso)]
            dfall_qsos[flux].plot(label='%s  , Chisq=%s'%(fqso.replace('.fits',''), dic_chisq[fqso]),
                                         xlim=xlimits, ylim=ylimits, ax=ax)
        plt.legend(loc='best')

        ax2 = plt.subplot(1, 2, 2)
        dfall_qsos[self.coadd_id].plot(label=self.coadd_id, xlim=xlimits, ylim=ylimits, ax=ax2)
        plt.legend(loc='best')
        plt.title('THING_ID: %s'%(self.th_id))
        plt.show(block=True)
        return 0



    def plot_chisq_dist(self, frac):
        """Save chisq for all spectra, and plot a histogram on the fly"""

        #for i in zipchisq.values():
        self.chisq_dist.append(frac)

        plt.hist(self.chisq_dist, bins=10, range=(0.0, 1.0))
        plt.ylabel('#')
        plt.xlabel('Fraction Chisq')
        plt.title('Chisq Histogram')
        plt.show(block=True)
        #plt.savefig('chisq.pdf')
        return 0



    def write_stats_open(self, rank):
        """Write all chisq and trim after eliminating trim_chisq > chisq"""
        self.write_stats = {'dist' : open(self.stats_file  + self.suffix.format(rank), 'w'),
                            'bad' : open(self.stats_file3 + self.suffix.format(rank), 'w')}
        self.write_stats['bad'].write('#Spec with pure noise: Chisq<2 \n')
	self.write_stats['bad'].write('#THING_ID, file, chisq \n')



    def write_stats_close(self):
        """Close files"""
        for i in ['dist', 'bad']:
            self.write_stats[i].close()



    def write_stats_file(self, frac, name):
        """Write chisq"""
        self.write_stats[name].write(str(frac) + '\n')
        self.write_stats[name].flush()



    def write_fits(self, result, lpix):
         data    = pd.concat([r for r in result], axis=1).fillna(0)
         nrows   = len(data.index)
         data    = data.reset_index().to_dict(orient='list')
         names   = list(data)
         formats = ['f4']*len(names)
         fdata   = np.zeros(nrows, dtype=dict(names= names, formats=formats))
         fdata   = {key:np.array(value) for key,value in data.items()}

         file_name = os.path.join(self.pix_dir, 'pix_%s.fits'%(lpix))
         if os.path.isfile(file_name): os.system("rm %s"%(file_name))
         fits = fitsio.FITS(file_name, 'rw')
         fits.write(fdata, header={'Healpix':'%s'%(lpix), self.full_file: 'Npix_side =%s'%(self.Npix_side) })
         fits.close()

         if self.verbose: print ('... Writing FITS file: %s'%(lpix))
	



    def extra_info(self, file, th_id, name):
        if file:
            _, plate, mjd, fiber = file.replace('.fits','').split('-')
            condition = 'THING_ID == %s & PLATE == %s & MJD == %s'%(th_id, plate, mjd)
            return float(self.df_qsos.query(condition)[name].values)
        else:
            return 0



    def resize_arr(self, arr):
	resize_str = np.array(arr)
	resize_str.resize(self.max_str)
	return resize_str



    def master_fits(self, lpix, thid, dict_z):
	#fits colums should have the same lenght :/
        self.max_str = max([len(i) for i in dict_z])
        all_lfiles=[]; z_war = []; z_err = []; z = []
	
        for dicts in dict_z:
            all_lfiles.append(self.resize_arr(dicts.keys()))
	    z_war.append(self.resize_arr([val[0] for val in dicts.values()]))
	    z_err.append(self.resize_arr([val[1] for val in dicts.values()]))
	    z.append(self.resize_arr([val[2] for val in dicts.values()]))
		
        data = {}
        data['healpix']   = np.array(lpix,            dtype='i4')
        data['thing_id']  = np.array(thid,            dtype='i4')
        data['files']     = np.array(all_lfiles,      dtype='S32')
        data['z_warning'] = np.array(z_war,           dtype='i4')
        data['z_err']     = np.array(z_err,           dtype='f4')
        data['z']         = np.array(z,               dtype='f4')

        os.system("rm master_file.fits")
        if self.verbose: print (' ... Writing MASTER-FITS file')
        fits = fitsio.FITS('master_file.fits','rw')
        fits.write(data)
        fits.close()


    def plot_stats(self, size):
        """Collect all chisq and plot a histogram"""

	chisq_ = 'chisq'
        total_chisq = []
        for i in np.arange(size):
            Chisq     = pd.read_csv(self.stats_file  + self.suffix.format(i), sep='\s', 
			names=['THING_ID','#number',chisq_] ,header=None)
            total_chisq.append(Chisq)
	
        df = pd.concat(total_chisq) 
 	tot_spec = df['#number'].sum()
        if self.verbose:
            print("\n Statistics of chisq distribution")
            print (df[chisq_].describe())
	
        bins  = 10
        range = (0,1)

        if self.use_bokeh:
            try:
                #testing different visualization methods
                from bokeh.plotting import figure, show

                TOOLS = "pan,box_zoom,reset,tap,save,crosshair"
                p1 = figure(title="Counting chisq for repeated THING_ID", tools=TOOLS)

                hist, edges = np.histogram(df[chisq_], density=False, bins=bins, range=range)
                p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
                        fill_color="blue", line_color="#FF7373", legend="Chisq-all: %s"%(len(df)))

                p1.xaxis.axis_label = 'Chisq'
                p1.yaxis.axis_label = '#'
                #output_file('histogram.html', title="histogram.py example")
                show(p1)
            except:
                print ('Install Bokeh for fun')
        else:
	    dict={}
	    for i in np.arange(1,19):
		x = np.array(np.histogram(df[df['#number']==i][chisq_].values, range=[0,1], bins=10)[0])
		dict[i] = np.array([np.log(j) if j!=0 else 0 for j in x ])
	    final = pd.DataFrame(index =[i/10. for i in np.arange(10)], data=dict)
	    print(final.sort(axis=1))
	
	    plt.figure(figsize = (18, 8))
	    ax = plt.subplot(111)
	    sns.heatmap(final.sort(axis=1), linewidths=0.5, annot=True, fmt=".1f", 
		linecolor='white', cmap="YlGnBu",  label='Total:%s'%(len(df)), ax=ax)
	    plt.ylabel('% accepted per THING_ID')	
	    plt.xlabel('Repeated THING_ID')  
	    plt.title('Total Spec:%s,   Unique THING_ID:%s'%(tot_spec, len(df)))
	    plt.legend(loc = 'best')
            plt.show(block=True)




if __name__=='__main__':
    print ("goofing around :P ")
    Qsos    = Qso_catalog(None, verbose = True)
    Qsos.plot_stats(8)

