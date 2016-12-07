
"""
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
    for i, lpix in enumerate(pixel[:]):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        if Qsos.verbose and i % 5 == 0: print ('#-- pix', i, {lpix: thingid_repeat})
        result = []

        for th_id in thingid_repeat.keys():
            dict_qso = Qsos.get_files(thing_id = th_id)
            old_qsos = len(list(dict_qso.keys()))
            dict_file, dict_chisq, dict_qso = Qsos.cal_chisq(dict_qso)

            len_files = len(list(dict_qso.keys()))
            if Qsos.write_hist:   Qsos.write_stats_file('dist', th_id, old_qsos, len_files)

            if len_files == 0: continue
            dfall_qsos = Qsos.coadds(dict_file)

            if Qsos.write_ffits:  result.append((th_id, dfall_qsos))
            if Qsos.show_plots:   Qsos.plot_coadds(dict_chisq)
            if Qsos.write_master:
                for info in dict_qso.values(): Qsos.all_info.append((lpix, th_id, info))

        if Qsos.write_ffits: Qsos.write_fits(result, lpix)
    return 0
