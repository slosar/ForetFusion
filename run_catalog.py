
"""
 using mpi4py

** TODO:
    - self.THING_ID, flux, ivar
    - add comments, and del sdss_catalog
    - README has some typos, and difficult to read
"""


import math
from qso_catalog import Qso_catalog
from main_file import split_pixel
from get_files import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dir_files = 'data/'
file_name = 'subset_spAll-v5_10_0.csv'

if rank == 0:
    df_fits = read_sub_fits(dir_files, file_name)
    Qsos    = Qso_catalog(df_fits)

    Qsos.rep_thid    = 4
    Qsos.verbose     = True
    Qsos.show_plots  = False
    Qsos.write_names = False
    Qsos.write_hist  = False

    Qsos.filtering_qsos(condition= Qsos.condition)
    unique_pixels = Qsos.adding_pixel_column()
    Qsos.ask_for_files(get_files= False)
    #print (Qsos.df_qsos.query('PIX == 6219 & (THING_ID == 77964771)'))
    if Qsos.write_names: Qsos.write_file_names()

    lenpix = len(unique_pixels)
    nchunk = int(math.ceil(lenpix*1./size))
    chunks = [unique_pixels[i:i+ nchunk] for i in range(0, lenpix, nchunk)]
else:
    chunks = []
    Qsos   = None


Qsos = comm.bcast(Qsos, root=0)
if Qsos.write_hist: Qsos.write_stats_open(rank)


chunk_pix  = comm.scatter(chunks, root=0)
split_pixel(chunk_pix, Qsos)
comm.Barrier()


if Qsos.write_hist: Qsos.write_stats_close()
if rank == 0 and Qsos.write_hist and Qsos.show_plots:
    Qsos.plot_stats(size)
    print ('-- stats are on Chisq_dist.csv files')
