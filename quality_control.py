import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook


summary = pd.read_csv('summary.txt', sep='\t')

#print len(phosphosites.index)
filtsum = summary[['Experiment','MS/MS Identified']]
#sumbrief2 = sumbrief.rows[0:50]
sumbrief = filtsum[0:50]

wcls = sumbrief[sumbrief['Experiment'].str.contains('_WCL')]
wcl = wcls[~wcls['Experiment'].str.contains('_WCLP')]
wclp = wcls[wcls['Experiment'].str.contains('_WCLP')]
ubs = sumbrief[sumbrief['Experiment'].str.contains('_Ub')]
ub = ubs[~ubs['Experiment'].str.contains('_UbP')]
ubp = ubs[ubs['Experiment'].str.contains('_UbP')]

#print wcl
#print wclp
#print ub
#print ubp


# compute the boxplot stats
ubstats = cbook.boxplot_stats(ub[["MS/MS Identified"]].values, whis='range', bootstrap=None, labels=None)
ubpstats = cbook.boxplot_stats(ubp[["MS/MS Identified"]].values, whis='range', bootstrap=None, labels=None)
wclstats = cbook.boxplot_stats(wcl[["MS/MS Identified"]].values, whis='range', bootstrap=None, labels=None)
wclpstats = cbook.boxplot_stats(wclp[["MS/MS Identified"]].values, whis='range', bootstrap=None, labels=None)

fs = 10 # fontsize

# demonstrate how to toggle the display of different elements:
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(4,4))
axes[0, 0].bxp(ubstats)
axes[0, 0].set_title('ub', fontsize=fs)

axes[0, 1].bxp(ubpstats)
axes[0, 1].set_title('ubp', fontsize=fs)

axes[1, 0].bxp(wclstats)
axes[1, 0].set_title('wcl', fontsize=fs)

axes[1, 1].bxp(wclpstats)
axes[1, 1].set_title('wclp', fontsize=fs)

data_x1 = [1, 1]
data_yub = [3615, 2430]
data_yubp = [2868, 2401]
data_ywcl = [7099, 7197]
data_ywclp = [9831, 6291]

axes[0, 0].plot(data_x1, data_yub, "or")
axes[0, 1].plot(data_x1, data_yubp, "or")
axes[1, 0].plot(data_x1, data_ywcl, "or")
axes[1, 1].plot(data_x1, data_ywclp, "or")

for ax in axes.flatten():
    ax.set_yticklabels(0, 10000, 500)

fig.subplots_adjust(hspace=0.4)
plt.show()