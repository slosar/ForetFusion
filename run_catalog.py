
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
    with open('Chisq_dist_{}.csv'.format(rank), 'w') as write_stats, \
    open('Chisq_dist_sec_{}.csv'.format(rank), 'w') as write_stats_two:
    print 'Number of healpix:', len(pixel)
        for i, lpix in enumerate(pixel):
            thingid_repeat = Qsos.pix_uniqueid(lpix)
            if not thingid_repeat: continue
            if i%2 ==0: print i, {lpix: thingid_repeat}

            for thids in thingid_repeat:
                write_stats.flush(), write_stats_two.flush()
                qso_files = Qsos.get_files(thing_id= thids)

                flag = 1
                while flag:
                    dfall_qsos = Qsos.coadds(qso_files, Qsos.spec_cols)
                    zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                    #some plots
                    #Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
                    #if flag==1: Qsos.plot_chisq_dist(zipchisq)

                    if flag==1:
                        for chi in zipchisq.values(): write_stats.write(str(chi) + '\n')
                 
                    #check specs that have chisq > self.del_chisq, if none, get out
                    flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Qsos.del_chisq))
                    if flag==0:
                        for chi in zipchisq.values(): write_stats_two.write(str(chi) + '\n')
                        continue
	
                    qso_files = Qsos.select_chisq(zipchisq, Qsos.del_chisq)
                    if len(qso_files) == 0:
                        print 'Really bad measurement, THING_ID:', thids
                        flag=0



Pars    = Ini_params()
df_fits = read_subset_fits(Pars.dir_fits, Pars.sub_file)
Qsos    = Qso_catalog(df_fits)
Qsos.filtering_qsos(condition= Pars.condition)
unique_pixels = Qsos.adding_pixel_column()

#test
#print Qsos.df_qsos.query('PIX == 6219 & (THING_ID == 77964771 | THING_ID== 68386221)')
Qsos.ask_for_files(get_them= False)

if rank == 0:
    #for targ, bits in Qsos.targets.iteritems():
    #    print "Quasars in {} =".format(targ), Qsos.searching_quasars(targ, bits).sum()
    #Qsos.print_file_names()

    lpix = len(unique_pixels)
    n = int(math.ceil(lpix*1./size))
    chunks = [unique_pixels[i:i+n] for i in range(0, lpix, n)]
else:
    chunks = []


chunk_pix = comm.scatter(chunks, root=0)
split_pixel(chunk_pix, Qsos, rank)
comm.Barrier()

print 'stats on Chisq_dist_.csv files' 
