import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
proteins = pd.read_table('proteinGroups.txt', low_memory=False)

intensity_cols = [c for c in proteins.columns if 'intensity '
            in c.lower() and 'lfq' not in c.lower()]


wcl_cols = [c for c in intensity_cols if '_wcl' in c.lower() and '_wclp' not in c.lower()]
wclp_cols = [c for c in intensity_cols if '_wclp' in c.lower()]
ub_cols = [c for c in intensity_cols if '_ub' in c.lower() and '_ubp' not in c.lower()]
ubp_cols = [c for c in intensity_cols if '_ubp' in c.lower()]

mask = (proteins['Reverse'] != '+') & \
       (proteins['Potential contaminant'] != '+')

intensities = proteins[mask][intensity_cols]
total_intensities = proteins[intensity_cols].sum(axis=0)
normed_intensities = intensities / total_intensities

idx = (normed_intensities != 0).any(axis=1)

names = proteins[mask][idx]['Protein IDs']
nonzero_intensities = normed_intensities[idx]

wcl = nonzero_intensities[wcl_cols]
wclp = nonzero_intensities[wclp_cols]
ub = nonzero_intensities[ub_cols]
ubp = nonzero_intensities[ubp_cols]

wcl_ctrl = [c for c in wcl.columns if 'control' in c.lower()]
wclp_ctrl = [c for c in wclp.columns if 'control' in c.lower()]
ub_ctrl = [c for c in ub.columns if 'control' in c.lower()]
ubp_ctrl = [c for c in ubp.columns if 'control' in c.lower()]

wcl_exp = [c for c in wcl.columns if 'control' not in c.lower()]
wclp_exp = [c for c in wclp.columns if 'control' not in c.lower()]
ub_exp = [c for c in ub.columns if 'control' not in c.lower()]
ubp_exp = [c for c in ubp.columns if 'control' not in c.lower()]

# Need to use underlying numpy arrays for singleton expansion ('broadcasting')
# and form new DataFrame using appropriate column names.
wcl_foldch = pd.DataFrame(np.log2(wcl[wcl_exp]).values - np.log2(wcl[wcl_ctrl]).values, columns=wcl_exp)
wclp_foldch = pd.DataFrame(np.log2(wclp[wclp_exp]).values - np.log2(wclp[wclp_ctrl]).values, columns=wclp_exp)
ub_foldch = pd.DataFrame(np.log2(ub[ub_exp]).values - np.log2(ub[ub_ctrl]).values, columns=ub_exp)
ubp_foldch = pd.DataFrame(np.log2(ubp[ubp_exp]).values - np.log2(ubp[ubp_ctrl]).values, columns=ubp_exp)


# 2nd-to-last element is Whangee TPK1-KO and penultimate element is Whangee Tunicamycin.
# Only histogram finite (non-inf, non-NaN) values.

tuni = [wcl_foldch[wcl_foldch.columns[-1]][np.isfinite(wcl_foldch[wcl_foldch.columns[-1]])].values, wclp_foldch[wclp_foldch.columns[-1]][np.isfinite(wclp_foldch[wclp_foldch.columns[-1]])].values, ub_foldch[ub_foldch.columns[-1]][np.isfinite(ub_foldch[ub_foldch.columns[-1]])].values, ubp_foldch[ubp_foldch.columns[-1]][np.isfinite(ubp_foldch[ubp_foldch.columns[-1]])].values]

plt.hist(tuni, normed=1, histtype='bar',
                            color=['b', 'r', 'y', 'g'],
                            label=['WCL', 'WCLP', 'UB', 'UBP'])

plt.legend()

plt.title('Tunicamycin')

plt.show()


tpk = [wcl_foldch[wcl_foldch.columns[-2]][np.isfinite(wcl_foldch[wcl_foldch.columns[-2]])].values, wclp_foldch[wclp_foldch.columns[-2]][np.isfinite(wclp_foldch[wclp_foldch.columns[-2]])].values, ub_foldch[ub_foldch.columns[-2]][np.isfinite(ub_foldch[ub_foldch.columns[-2]])].values, ubp_foldch[ubp_foldch.columns[-2]][np.isfinite(ubp_foldch[ubp_foldch.columns[-2]])].values]

plt.hist(tpk, normed=1, histtype='bar',
                            color=['b', 'r', 'y', 'g'],
                            label=['WCL', 'WCLP', 'UB', 'UBP'])

plt.legend()

plt.title('TPK1 KO')
plt.show()

#Adds protein IDs to the dataframes
ubp_foldch['names'] = names
ubp_foldch.set_index(names, inplace=True)

ub_foldch['names'] = names
ub_foldch.set_index(names, inplace=True)

wcl_foldch['names'] = names
wcl_foldch.set_index(names, inplace=True)

wclp_foldch['names'] = names
wclp_foldch.set_index(names, inplace=True)

#writes CSV files with Protein IDs sorted by fold change
idx = ubp_foldch[ubp_foldch.columns[-2]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./ubp_tuni_proteinID_fold-change_rank.csv', sep = '\t')

idx = ub_foldch[ub_foldch.columns[-2]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./ub_tuni_proteinID_fold-change_rank.csv', sep = '\t')


idx = wcl_foldch[wcl_foldch.columns[-2]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./wcl_tuni_proteinID_fold-change_rank.csv', sep = '\t')


idx = wclp_foldch[wclp_foldch.columns[-2]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./wclp_tuni_proteinID_fold-change_rank.csv', sep = '\t')

idx = ubp_foldch[ubp_foldch.columns[-3]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./ubp_tpk1_proteinID_fold-change_rank.csv', sep = '\t')

idx = ub_foldch[ub_foldch.columns[-3]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./ub_tpk1_proteinID_fold-change_rank.csv', sep = '\t')


idx = wcl_foldch[wcl_foldch.columns[-3]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./wcl_tpk1_proteinID_fold-change_rank.csv', sep = '\t')

idx = wclp_foldch[wclp_foldch.columns[-3]]
with pd.option_context('mode.use_inf_as_null', True):
    idx = idx.dropna()
obj = idx[idx.argsort()]
obj.to_csv('./wclp_tpk1_proteinID_fold-change_rank.csv', sep = '\t')
