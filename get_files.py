#!/usr/bin/env python

import os
import fitsio
import pandas as pd
from mpi4py import MPI


def read_sub_fits(dir_fits, file_name):
    """Read the subsample of spAll we are interested in. Return a DataFrame"""
    file_name = os.path.join(dir_fits, file_name)
    if not os.path.isfile(file_name):
        print ('File not found: {}'.format(file_name))
        MPI.COMM_WORLD.Abort()

    print ('Reading file...')
    df = pd.read_csv(file_name, sep=',')
    return df



def read_fits(dir_fits, file_name, fits_cols):
    """Read selected columns in the .fits file. Return a DataFrame"""
    file_name = os.path.join(dir_fits, file_name)
    if not os.path.isfile(file_name):
        print('File not found: {}'.format(file_name))
        MPI.COMM_WORLD.Abort()

    fits       = fitsio.FITS(file_name)

    #http://stackoverflow.com/questions/30283836/
    # creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors
    fits_read  = fits[1].read(columns= fits_cols)
    fits_to_df = {col:fits_read[col].byteswap().newbyteorder() for col in fits_cols}
    df = pd.DataFrame(fits_to_df)
    return df



if False:
    #Just In case we need to get the spec files from the sdss-website
    #we need import mechanize and from base64 import b64encode
    def get_web_files(self, file_name, passwd):
        """nasty hack to get files from sdss website, but will change it later,
            only works with Python2"""

        self.sdss_url    = 'https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spectra/'
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


