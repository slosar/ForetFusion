
"""
Same as sdss_catalog.py, but no comments

Running the qso_catalog.py:
* It goes through the SpAll file and filter all qsos in the BOSS, EBOSS
that satisfy the bit condition and also
    'CLASS== "QSO" & (OBJTYPE=="QSO" | ''OBJTYPE=="NA".) & THING_ID != -1'
* Computes the healpix given the RA and DEC, and group all of them by healpix.
* Withing a healpix, finds the objects with the same THING_ID
* Coadd them (average ivar*flux.) and compute the chisq
* Loop over:  If the chisq is more than 4, eliminate that spec,
    coadd again and get new chisq

"""

def split_pixel(pixel, Qsos):
    for i, lpix in enumerate(pixel[:10]):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        if Qsos.verbose and i % 5 == 0: print (i, {lpix: thingid_repeat})

        for thids in thingid_repeat:
            qso_files = Qsos.get_files(thing_id = thids)

            flag = 1
            while flag:
                dfall_qsos = Qsos.coadds(qso_files)
                zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                if Qsos.write_stats and flag == 1:
                    Qsos.write_stats_file(zipchisq, 'all')
                    Qsos.write_stats['all'].flush(), Qsos.write_stats['trim'].flush()

                if Qsos.show_plots:
                    Qsos.plot_coadds(dfall_qsos, thids, zipchisq)
                    if flag == 1: Qsos.plot_chisq_dist(zipchisq)


                #check specs that have chisq > self.trim_chisq, if none, get out
                flag = len(qso_files) - len(Qsos.ftrim_chisq(zipchisq))
                if flag == 0:
                    if Qsos.write_stats: Qsos.write_stats_file(zipchisq, 'trim')
                    continue

                qso_files = Qsos.ftrim_chisq(zipchisq)
                if len(qso_files) == 0:
                    flag=0
                    print ('Really bad measurement, THING_ID:', thids)
