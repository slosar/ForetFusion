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


![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/THING_ID_4.jpg)



* Write a file for each healpix that contains THING_ID's. For each THING_ID coadd 'flux'. For example, considering the Healpix 24985, there are three different THING_ID with their corresponding coadd, and_mask and or_mask:
	
![](https://github.com/ja-vazquez/ForetFusion/blob/master/figs/healpix.jpg)









