
"""
Same as run_catalog.py, but no comments and using mpi4py
"""

import math
from qso_catalog import *
from get_files import *
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def split_pixel(pixel, Qsos, rank):
    write_files   =  open('SpAll_files_{}.csv'.format(rank),'w')
    write_stats   =  open('Chisq_dist_{}.csv'.format(rank), 'w')
    write_stats_2 =  open('Chisq_dist_sec_{}.csv'.format(rank), 'w')

    for i, lpix in enumerate(pixel):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue

        if i%10 ==0: print {lpix: thingid_repeat}
        write_files.flush(); write_stats.flush(), write_stats_2.flush()

        for thids in thingid_repeat:
            qso_files = Qsos.get_files(thing_id= thids)
            for names in qso_files: write_files.write('v5_10_0/spectra/' + names + '\n')

            flag, flag_2 = 1, 0
            while flag:
                dfall_qsos = Qsos.coadds(qso_files, Qsos.spec_cols)
                zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                #some plots
                #Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
                #if flag==1: Qsos.plot_chisq_dist(zipchisq)

                if flag_2:
                    for chi in zipchisq.values(): write_stats_2.write(str(chi) + '\n')
                if flag:
                    for chi in zipchisq.values(): write_stats.write(str(chi) + '\n')
                    flag_2 +=1


                #check specs that have chisq > self.del_chisq, if none, get out
                flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Qsos.del_chisq))
                qso_files = Qsos.select_chisq(zipchisq, Qsos.del_chisq)

    write_files.close()
    write_stats.close()
    write_stats_2.close()




Pars    = Ini_params()
df_fits = read_subset_fits(Pars.dir_fits, Pars.sub_file)
Qsos    = Qso_catalog(df_fits)
Qsos.filtering_qsos(condition= Pars.condition)
unique_pixels = Qsos.adding_pixel_column() #[:100]
Qsos.need_files = '1'


if rank == 0:
    lpix = len(unique_pixels)
    n = int(math.ceil(lpix*1./size))
    chunks = [unique_pixels[i:i+n] for i in range(0, lpix, n)]
else:
    chunks = []


chunk_pix = comm.scatter(chunks, root=0)
split_pixel(chunk_pix, Qsos, rank)
comm.Barrier()

