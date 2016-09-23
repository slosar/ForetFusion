
"""
 using mpi4py

** TODO:
    - self.THING_ID, flux, ivar
    - add comments, and del sdss_catalog
    - README has some typos, and difficult to read
"""


from qso_catalog import Qso_catalog
import matplotlib.pyplot as plt
from get_files import *
import numpy as np

split     = 20
dir_files = 'data/'
file_name = 'subset_spAll-v5_10_0.csv'


df_fits = read_sub_fits(dir_files, file_name)
Qsos    = Qso_catalog(df_fits, verbose = False)

Qsos.my_own_filter(condition= 'OBJTYPE=="SKY".ljust(16)')
qso_files = Qsos.write_file_names()
qso_test  = qso_files.tolist()[:400]
[(Qsos.get_bnl_files(qso.split('/')[0],qso) if not os.path.isfile(Qsos.dir_spec + qso) else print('have it')) for qso in qso_test]
qso_split = np.split(np.array(qso_test), split)

stack_flux     = []; stack_ivar     = []
for i, qso in enumerate(qso_split):
    Qsos.th_id = i
    dfall_qsos = Qsos.coadds(qso_split[i])
    stack_flux.append(dfall_qsos[Qsos.flux_ivar_id])
    stack_ivar.append(dfall_qsos[Qsos.ivar_id])

    zipchisq   = Qsos.calc_chisq(qso, dfall_qsos)
    Qsos.plot_coadds(dfall_qsos, zipchisq)

df_flux = (pd.concat([stack_flux[i] for i, _ in enumerate(qso_split)], axis=1).fillna(0))
df_ivar = (pd.concat([stack_ivar[i] for i, _ in enumerate(qso_split)], axis=1).fillna(0))
final   = pd.concat([df_flux.cumsum(axis=1).iloc[:,-1], df_ivar.cumsum(axis=1).iloc[:,-1]], axis=1)

final['coadd'] = final['flux_ivar_%s'%(split-1)]/final['ivar_%s'%(split-1)]

plt.figure(figsize = (18, 8))
final['coadd'].plot(xlim=[3.55, 4], ylim=[-10, 25], label= '%s files used'%(len(qso_test)))
plt.legend(loc='best')
plt.show(block=True)




