
import os, sys
import fitsio
import pandas as pd


def read_sub_fits(dir_fits, sub_file):
    """Read the subsample of spAll we are interested in. Return a DataFrame"""
    file_name = '{}{}'.format(dir_fits, sub_file)
    if not os.path.isfile(file_name):
        sys.exit('File amm not found: {}'.format(file_name))

    print ('Reading file...')
    df = pd.read_csv(file_name, sep=',')
    return df



def read_fits(dir_fits, file_name, columns):
    """Read selected columns in the spAll file. Return a DataFrame"""
    file_name = '{}{}.fits'.format(dir_fits, file_name)
    #if not os.path.isfile(file_name):
    #    sys.exit('File mmm not found: {}'.format(file_name))
     
    fits          = fitsio.FITS(file_name)
    fits_columns  = columns

    #http://stackoverflow.com/questions/30283836/
    # creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors
    fits_read  = fits[1].read().byteswap().newbyteorder()
    fits_to_df = {col: fits_read[col] for col in fits_columns}
    df = pd.DataFrame(fits_to_df)
    return df






