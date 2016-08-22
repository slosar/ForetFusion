# ForetFusion
Code to compute coadds for Lya-forest data in BOSS and eBOSS. 

Code:

    > mpirun -np 4 run_catalog.py

The calculations and main functions are located at:

     qso_catalog.py
     
Python libraries you may need: fitsio, pandas, mpi4py, healpy. I am using Python 3.5.2.

* It goes through the [SpAll.fits](https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spAll-v5_10_0.fits) file and filter all qsos in the BOSS, EBOSS (The filtering conditions can be easily modified).
(the full 12GB file can be read, however to save some time we are using a subset, [subset_spAll](http://www.cosmo.bnl.gov/www/jvazquez/try_1/), with the columns we are interested on)

* Filter the Quasars that satisfy the the bit condition: 

    > ('Quasars with Bit condition in EBOSS_TARGET1 =', 201742)
    
    > ('Quasars with Bit condition in EBOSS_TARGET0 =', 54108)
    
    > ('Quasars with Bit condition in BOSS_TARGET1 =' , 466098)
    
    > Total qsos = 3008000

* Filter the Quasars that satisfy the condition  'CLASS== "QSO" & (OBJTYPE=="QSO" | ''OBJTYPE=="NA".) & THING_ID != -1'
    
    > ('Quasars with General condition in EBOSS_TARGET1 =', 160854)
    
    > ('Quasars with General condition in EBOSS_TARGET0 =', 38164)
    
    > ('Quasars with General condition in BOSS_TARGET1 =' , 285071)
    
    > Total qsos = 610184
    
    
* Filtering by using both conditions, we have:

    > ('Quasars with Both condition in EBOSS_TARGET1 =', 160854)
    
    > ('Quasars with Both condition in EBOSS_TARGET0 =', 38164)
    
    > ('Quasars with Both condition in BOSS_TARGET1 = ' , 285071)
    
    > Total qsos = 484089

* The names of the files that contain the 484089 spectra corresponding to these quasars are located at [SpAll_files](https://github.com/ja-vazquez/ForetFusion/blob/master/SpAll_files.csv).
    If you want to write a file with all quasars names, or with some other condition: use 
    **Qsos.write_names = True**

* Once filtered and selected the Qsos, it computes the healpix given RA and DEC, and group quasar observations by healpix.
i.e. 

![](https://github.com/ja-vazquez/ForetFusion/blob/master/THING_ID_3.jpg )


* Withing a healpix, finds the repeated objects with the same THING_ID (can be changed in **self.rep_thid = 4**), Coadd them (average ivar*flux.) and compute the chisq:

N.B. In case you don't have the spec files, there are two options, (bnl): get them from the bnl cluster or (sdss): directly from the sdss website (although is password protected), by modifying **get_files= True**
i.e. 

![](https://github.com/ja-vazquez/ForetFusion/blob/master/THING_ID_1.jpg)

* Loop over:  If the chisq is more than 4 (or change this number in **self.trim_chisq = 4**), eliminate that spec,
    coadd again and get new chisq

![](https://github.com/ja-vazquez/ForetFusion/blob/master/THING_ID_2.jpg)


* Write a file for each healpix that contains THING_ID's and for each one the coadd 'flux','loglam':
	
		HEALPIX: THING_ID_1	:	'flux_1','loglam_1'  .fits/.csv
				 THING_ID_2	:	'flux_2','loglam_2'  .fits/.csv
				 THING_ID_3	:	'flux_3','loglam_3'  .fits/.csv
				 	.

* It can show some coadd and/or chis-distribution plots on the fly **Qsos.show_plots = True**.
* At the end shows a distribution of chisq taking into account all the Quasars and after eliminating
those with chisq >= 4. 

The plot below shows the distribution of chisq considering only the quasars with repeated THING_ID, in this case for testing, with more than 4 times.

![](https://github.com/ja-vazquez/ForetFusion/blob/master/chisq_4.jpg )



