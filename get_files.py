#!/usr/bin/env python

# Add concepto: deposito, retiro, iva

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



def read_fits(dir_fits, file_name, columns):
    """Read selected columns in the .fits file. Return a DataFrame"""
    file_name = os.path.join(dir_fits, file_name)
    if not os.path.isfile(file_name):
        print('File not found: {}'.format(file_name))
        MPI.COMM_WORLD.Abort()

    fits          = fitsio.FITS(file_name)
    fits_columns  = columns

    #http://stackoverflow.com/questions/30283836/
    # creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors
    fits_read  = fits[1].read().byteswap().newbyteorder()
    fits_to_df = {col: fits_read[col] for col in fits_columns}
    df = pd.DataFrame(fits_to_df)
    return df






