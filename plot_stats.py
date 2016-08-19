
import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.mpl_style', 'default')
#matplotlib.style.use('ggplot')


file_name = 'Chisq_dist'
num_files = 4

total_chisq = []
total_chisq_sec = []
for i in range(num_files):
    Chisq = pd.read_csv(file_name + '_{}.csv'.format(i))
    Chisq_sec = pd.read_csv(file_name + '_sec_{}.csv'.format(i))
    total_chisq.extend(Chisq.values.flatten())
    total_chisq_sec.extend(Chisq_sec.values.flatten())

plt.figure()
ax = plt.subplot(111)
df = pd.DataFrame(total_chisq, columns=['chisq'])
df_sec = pd.DataFrame(total_chisq_sec, columns=['chisq'])

df['chisq'].plot.hist(bins=80, range=(0.1, 6), alpha=0.9, ax=ax, color='r', label='chisq, %s'%(len(df)))
df_sec['chisq'].plot.hist(bins=80, range=(0.1, 6), alpha=0.5, ax=ax, color='b', label='after loop, %s'%(len(df_sec)))

plt.ylabel('#')
plt.xlabel('chisq')
plt.legend(loc = 'best')

print (df['chisq'].describe())
print (df_sec['chisq'].describe())

plt.show(block=True)
