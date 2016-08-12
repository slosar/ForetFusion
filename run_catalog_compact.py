




def split_pixel(pixel, Qsos, write_var, rank):
    for _, lpix in enumerate(pixel):
        thingid_repeat = Qsos.pix_uniqueid(lpix)
        if not thingid_repeat: continue
        print 'pixel' ,lpix
        for thids in thingid_repeat:
            qso_files = Qsos.get_files(thing_id= thids)
            for names in qso_files: write_var[rank].write('v5_10_0/spectra/' + names + '\n')
            #flag = 1
            #while flag:
            #    dfall_qsos = Qsos.coadds(qso_files, Qsos.spec_cols)
            #    zipchisq   = Qsos.calc_chisq(qso_files, dfall_qsos)

                #check specs that have chisq > self.del_chisq, if none, get out
            #    flag = len(qso_files) - len(Qsos.select_chisq(zipchisq, Qsos.del_chisq))
            #    qso_files = Qsos.select_chisq(zipchisq, Qsos.del_chisq)



