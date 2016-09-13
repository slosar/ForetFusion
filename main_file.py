
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
    for i, lpix in enumerate(pixel):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        if Qsos.verbose and i % 5 == 0: print ('#pix', i, {lpix: thingid_repeat})

        result = []
        all_qso_files = {}
        for th_id in thingid_repeat:
            qso_files = Qsos.get_files(thing_id = th_id)

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
                    result.append(dfall_qsos[[Qsos.coadd_id, Qsos.ivar_id, Qsos.and_mask_id, Qsos.or_mask_id]])
                    thih = "th_%s"%(th_id)
                    all_qso_files[thih]= qso_files
                    if Qsos.write_hist:
                        Qsos.write_stats_file(zipchisq, 'trim')
                        Qsos.write_stats['all'].flush(), Qsos.write_stats['trim'].flush()
                    continue


                qso_files = Qsos.ftrim_chisq(zipchisq)
                if len(qso_files) == 0:
                    if Qsos.write_hist: Qsos.write_stats['bad'].write(str(Qsos.th_id) + '\n')
                    if Qsos.verbose: print ('Really bad measurement, THING_ID:', Qsos.th_id)
                    flag = 0

        if Qsos.write_hist:
            x = '{"%s":%s}\n'%(lpix,all_qso_files)
            Qsos.write_stats['master'].write(str(x).replace('\'','"'))
        Qsos.write_fits(result, all_qso_files, lpix)
