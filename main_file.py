
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
import pandas as pd

def split_pixel(pixel, Qsos):
    for i, lpix in enumerate(pixel[:1]):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        if Qsos.verbose and i % 5 == 0: print (i, {lpix: thingid_repeat})

        result= []
        for thids in thingid_repeat:
            qso_files = Qsos.get_files(thing_id = thids)

            flag = 99
            while flag:
                dfall_qsos = Qsos.coadds(qso_files)
                zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                #write stats files and show some plots
                if flag == 99:
                    if Qsos.write_hist: Qsos.write_stats_file(zipchisq, 'all')
                    if Qsos.show_plots: Qsos.plot_chisq_dist(zipchisq)
                if Qsos.show_plots: Qsos.plot_coadds(dfall_qsos, zipchisq)


                #check specs that have chisq < self.trim_chisq, if none, get out
                flag = len(qso_files) - len(Qsos.ftrim_chisq(zipchisq))

                if flag == 0:
                    result.append(dfall_qsos[[Qsos.coadd_id, Qsos.ivar_id]])
                    if Qsos.write_hist:
                        Qsos.write_stats_file(zipchisq, 'trim')
                        Qsos.write_stats['all'].flush(), Qsos.write_stats['trim'].flush()
                    continue


                qso_files = Qsos.ftrim_chisq(zipchisq)
                if len(qso_files) == 0:
                    flag = 0
                    Qsos.write_stats['bad'].write(str(Qsos.thids) + '\n')
                    print ('Really bad measurement, THING_ID:', Qsos.thids)
        print (pd.concat([r for r in result], axis=1).fillna(0).head(20))