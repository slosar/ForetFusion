# ForetFusion
Code to compute coadds for Lya-forest data in BOSS and eBOSS. 

Code:

    > mpirun -np 4 run_catalog.py

The calculations and main functions are located at:

     qso_catalog.py
     
Python libraries you may need: fitsio, pandas, mpi4py, healpy. I am using Python 3.5.2.

* It goes through the [spAll-v5_10_0.fits](https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_0/spAll-v5_10_0.fits) file and filter all qsos within BOSS, EBOSS (filtering conditions can be easily modified).
(the full 12GB file can be read, however to save some time we are using a subset, [subset_spAll](www.cosmo.bnl.gov/www/jvazquez/forest/subset_spAll-v5_10_0.csv), with only the columns we are interested on)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/Filters.jpg )

* Filter the Quasars that satisfy the the [bit condition](http://www.sdss.org/dr12/algorithms/bitmasks/#BOSSTILE_STATUS): 

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

* The names of the files that contain the 484089 spectra corresponding to these quasars are located at [SpAll_files.csv](www.cosmo.bnl.gov/www/jvazquez/forest/SpAll_files.csv), use 
    **Qsos.write_names = True**

* Once filtered and selected the Qsos, we compute the healpix given RA and DEC. Then group quasar observations given by healpix.
i.e. 

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID.jpg )


* Withing a healpix, finds the repeated objects with the same THING_ID (i.e. **self.rep_thid = 4**), Coadd them (average ivar*flux.) and compute the chisq with respect to zero flux:

* Eliminate those observations with Chisq< 2 (**trim_chisq  = 2**) to neglect noise, i.e.






![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_1info.jpg)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_1.jpg)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_2.jpg)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_3.jpg)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_4.jpg)



* Write a file for each healpix that contains THING_ID's and for each one the coadd 'flux','loglam':
	
		HEALPIX: THING_ID_1	:	'flux_1','loglam_1'  .fits/.csv
				 THING_ID_2	:	'flux_2','loglam_2'  .fits/.csv
				 THING_ID_3	:	'flux_3','loglam_3'  .fits/.csv
				 	.

* It can show some coadd and/or chis-distribution plots on the fly **Qsos.show_plots = True**.
* At the end shows a distribution of chisq taking into account all the Quasars and after eliminating
those with chisq >= 4. 

The plot below shows the distribution of chisq considering only the quasars with repeated THING_ID, in this case for testing, with more than 4 times.

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/chisq_4.jpg )


* Last, there are some very bad observations where most of the specs have a chisq>4 compared to the coadd (and I'm not so sure what to do in these cases)

![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_4.jpg )
