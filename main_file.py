
"""
Same as run_catalog.py, but no comments and
"""

def split_pixel(pixel, Qsos):
    for i, lpix in enumerate(pixel):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        if i % 5 == 0: print (i, {lpix: thingid_repeat})

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
                flag = len(qso_files) - len(Qsos.select_chisq(zipchisq))
                if flag == 0:
                    if Qsos.write_stats: Qsos.write_stats_file(zipchisq, 'trim')
                    continue

                qso_files = Qsos.select_chisq(zipchisq)
                if len(qso_files) == 0:
                    flag=0
                    print ('Really bad measurement, THING_ID:', thids)
