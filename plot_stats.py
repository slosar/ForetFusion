
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

total_chisq = []
total_chisq_sec = []
for i in range(4):
    Chisq = pd.read_csv('Chisq_dist_{}.csv'.format(i))
    Chisq_sec = pd.read_csv('Chisq_dist_sec_{}.csv'.format(i))
    total_chisq.extend(Chisq.values.flatten())
    total_chisq_sec.extend(Chisq_sec.values.flatten())

plt.figure()
ax = plt.subplot(111)
df = pd.DataFrame(total_chisq, columns=['chisq'])
df_sec = pd.DataFrame(total_chisq_sec, columns=['chisq'])

df['chisq'].plot.hist(bins=50, range=(0.1,6), alpha=0.5, ax=ax)
df_sec['chisq'].plot.hist(bins=50, range=(0.1,6), alpha=0.5, ax=ax)

print df['chisq'].describe()
print df_sec['chisq'].describe()

plt.show()