
"""
 using mpi4py
"""

import math
from qso_catalog import Qso_catalog
from mpi4py import MPI
from main_file import split_pixel
from get_files import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


dir_files = 'data/'
file_name = 'subset_spAll-v5_10_0.csv'

df_fits = read_sub_fits(dir_files, file_name)
Qsos    = Qso_catalog(df_fits)

Qsos.show_plots  = True
Qsos.write_names = False

Qsos.filtering_qsos(condition= Qsos.condition)
unique_pixels = Qsos.adding_pixel_column()
Qsos.ask_for_files(get_files= False)

#test
#print Qsos.df_qsos.query('PIX == 6219 & (THING_ID == 77964771 | THING_ID== 68386221)')

if rank == 0:
    if Qsos.write_names: Qsos.write_file_names()
    for targ, bits in Qsos.targets.items():
        print ("Quasars in {} =".format(targ), Qsos.searching_quasars(targ, bits).sum())

    lenpix = len(unique_pixels)
    nchunk = int(math.ceil(lenpix*1./size))
    chunks = [unique_pixels[i:i+ nchunk] for i in range(0, lenpix, nchunk)]
else:
    chunks = []


Qsos.write_stats_open(rank)

chunk_pix  = comm.scatter(chunks, root=0)
split_pixel(chunk_pix, Qsos)
comm.Barrier()

Qsos.write_stats_close()
print ('stats are on Chisq_dist.csv files')
