
import math
from qso_catalog import *
from get_files import *
from run_catalog_compact import *
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



Pars = Ini_params()
df_fits = read_subset_fits(Pars.dir_fits, Pars.sub_file)
Qsos    = Qso_catalog(df_fits)
Qsos.filtering_qsos(condition= Pars.condition)
unique_pixels = Qsos.adding_pixel_column()[:100]
Qsos.read_ask = '3'
write_var = [open('SpAll_files_{}.csv'.format(i),'w') for i in range(size)]
#partial_pix   = functools.partial(split_pixel, Qsos=Qsos)


if rank == 0:
    lpix = len(unique_pixels)
    n = int(math.ceil(lpix*1./size))
    chunks = [unique_pixels[i:i+n] for i in range(0, lpix, n)]
else:
    chunks = []


y2 = comm.scatter(chunks, root=0)
split_pixel(y2, Qsos, write_var, rank)
comm.Barrier()

