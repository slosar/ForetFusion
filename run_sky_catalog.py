
"""
** TODO:
    - self.THING_ID, flux, ivar
    - add comments, and del sdss_catalog
    - parallelize it
"""


from qso_catalog import Qso_catalog
import matplotlib.pyplot as plt
from get_files import *
import numpy as np
import time

split     = 200
num_spec  = 10000
dir_files = 'data/'
file_name = 'subset_spAll-v5_10_0.csv'


df_fits = read_sub_fits(dir_files, file_name)
Qsos    = Qso_catalog(df_fits, verbose = False)

Qsos.my_own_filter(condition= 'OBJTYPE=="SKY".ljust(16)')
qso_files = Qsos.write_file_names()
qso_test = qso_files.tolist()[:num_spec]

#in case we need to get the files from the bnl, nasty but works
#[(Qsos.get_bnl_files(qso.split('/')[0],qso) if not os.path.isfile(Qsos.dir_spec + qso) else print('have it')) for qso in qso_test]
qso_split = np.split(np.array(qso_test), split)

start = time.time()
#here is where we parallelize
stack_flux     = []; stack_ivar     = []
for i, qso in enumerate(qso_split):
    Qsos.th_id = i
    dfall_qsos = Qsos.coadds(qso_split[i])
    stack_flux.append(dfall_qsos[Qsos.flux_ivar_id])
    stack_ivar.append(dfall_qsos[Qsos.ivar_id])

    #zipchisq   = Qsos.calc_chisq(qso, dfall_qsos)
    #Qsos.plot_coadds(dfall_qsos, zipchisq)

    if i%10==0: 
       print ('files', (i+1)*num_spec/split, 'timing', time.time() - start)


df_flux = (pd.concat([stack_flux[i] for i, _ in enumerate(qso_split)], axis=1).fillna(0))
df_ivar = (pd.concat([stack_ivar[i] for i, _ in enumerate(qso_split)], axis=1).fillna(0))
final   = pd.concat([df_flux.cumsum(axis=1).iloc[:,-1], df_ivar.cumsum(axis=1).iloc[:,-1]], axis=1)

final['coadd'] = final['flux_ivar_%s'%(split-1)]/final['ivar_%s'%(split-1)]

plt.figure(figsize = (18, 8))
final['coadd'].plot(ylim=[-10,100], label= '%s files used'%(len(qso_test)))
plt.legend(loc='best')
plt.show(block=True)


