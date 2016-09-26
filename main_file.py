
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
        if Qsos.verbose and i % 5 == 0: print ('#-- pix', i, {lpix: thingid_repeat})
        result = []

        for th_id in thingid_repeat.keys():
            qso_files = Qsos.get_files(thing_id = th_id)
            dic_file, dic_chisq = Qsos.cal_chisq(qso_files)

            frac = len(dic_file)/len(qso_files)

            if Qsos.write_hist: Qsos.write_stats_file(frac, 'dist')
            if Qsos.show_plots: Qsos.plot_chisq_dist(frac)

            if(frac == 0):
                if Qsos.write_hist: Qsos.write_stats['bad'].write(str(Qsos.th_id) + '\n')
                if Qsos.verbose: print ('Really bad measurement, THING_ID:', th_id)
                continue

            dfall_qsos = Qsos.coadds(dic_file)

            if Qsos.write_ffits:  result.append(dfall_qsos[[Qsos.coadd_id, Qsos.ivar_id, Qsos.and_mask_id,Qsos.or_mask_id]])
            if Qsos.write_master: Qsos.all_lpix.append(lpix); Qsos.all_thid.append(th_id);Qsos.all_qfiles.append(qso_files)
            if Qsos.show_plots: Qsos.plot_coadds(dfall_qsos, dic_chisq)

        if Qsos.write_ffits: Qsos.write_fits(result, lpix)
