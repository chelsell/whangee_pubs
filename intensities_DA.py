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

# 2nd-to-last element is Shmoo / CaCl2.
# Only histogram finite (non-inf, non-NaN) values.

wcltuni = plt.hist(wcl_foldch[wcl_foldch.columns[-1]][np.isfinite(wcl_foldch[wcl_foldch.columns[-1]])].values, label='WCL')

wclptuni = plt.hist(wclp_foldch[wclp_foldch.columns[-1]][np.isfinite(wclp_foldch[wclp_foldch.columns[-1]])].values, label='WCLP')

ubtuni = plt.hist(ub_foldch[ub_foldch.columns[-1]][np.isfinite(ub_foldch[ub_foldch.columns[-1]])].values, label='UB')

ubptuni = plt.hist(ubp_foldch[ubp_foldch.columns[-1]][np.isfinite(ubp_foldch[ubp_foldch.columns[-1]])].values, label='UBP')

plt.title('Tunicamycin')
plt.show()



wcltpk = plt.hist(wcl_foldch[wcl_foldch.columns[-2]][np.isfinite(wcl_foldch[wcl_foldch.columns[-2]])].values, label='WCL')

wclptpk = plt.hist(wclp_foldch[wclp_foldch.columns[-2]][np.isfinite(wclp_foldch[wclp_foldch.columns[-2]])].values, label='WCLP')

ubtpk = plt.hist(ub_foldch[ub_foldch.columns[-2]][np.isfinite(ub_foldch[ub_foldch.columns[-2]])].values, label='UB')

ubptpk = plt.hist(ubp_foldch[ubp_foldch.columns[-2]][np.isfinite(ubp_foldch[ubp_foldch.columns[-2]])].values, label='UBP')

plt.title('TPK1 KO')
plt.show()
