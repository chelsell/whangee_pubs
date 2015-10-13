import pandas
import cPickle as pic
import re
import numpy as np
import matplotlib.pyplot as plt
import urllib2
import collections

evidenceDataFrame = pandas.DataFrame.from_csv('evidence.txt', sep='\t')
proteinGroupsDataFrame = pandas.DataFrame.from_csv('proteinGroups.txt', sep='\t')
PhosphoSTYsitesDataFrame = pandas.DataFrame.from_csv('Phospho(STY)Sites.txt', sep='\t')
peptideDataFrame = pandas.DataFrame.from_csv('peptides.txt', sep='\t')

def updateDF(DataFrame, stringKey, stringValue):
    DataFrame = DataFrame[DataFrame[stringKey] != stringValue]

def denormalize(dfin, col):
    df = dfin
    drp = []
    new_entries = []
    for i in reversed(xrange(0, df.shape[0])):
        ids = df.iloc[i][col]
        if type(ids) is not str:
            continue
        tok = ids.split(";")
        if len(tok) > 1:
            idx = [] # needs indices of other columns with multiple values
            # for t in tok:
            for j in xrange(0, len(tok)):
                new_entry = df.iloc[i]
                new_entry[col] = tok[j]
                # for k in idx:
                    # set new_entry columns named in idx to df.iloc[i][k].split(";")[j]
                new_entries.append(new_entry)
            drp.append(i)
    df = df.drop(df.index[drp])
    df = df.append(new_entries)
    return df

def phosphoParser():
	phospho = pandas.read_csv("Phospho (STY)Sites.txt", sep = '\t', low_memory = False)
	deleteCol = []
	for i in range(len(phospho.columns)):
		if "UbP" in phospho.columns[i]:
			deleteCol.append(i)
	phospho = phospho.drop(phospho.columns[deleteCol], axis = 1)
	phospho = phospho[phospho['Reverse'] != '+']
	phospho = phospho[phospho['Potential contaminant'] != '+']
	return phospho

def proteinParser():
	protein = pandas.read_csv("proteinGroups.txt", sep = '\t', low_memory = False)
	deleteCol = []
	for i in range(len(protein.columns)):
		if "UbP" in protein.columns[i]:
			deleteCol.append(i)
	protein = protein.drop(protein.columns[deleteCol], axis = 1)
	protein = protein[protein['Reverse'] != '+']
	protein = protein[protein['Potential contaminant'] != '+']
	return protein

def intensity():
	file = proteinParser()
	file_unNorm = proteinParser()
	sum_col = file.sum(0)
	numIntensities = 0
	colIntensities = []
	ctrlIntensities = []
	for i in range(len(file.columns)):
		if "Intensity" in file.columns[i]:
			file[file.columns[i]] /= sum_col[i]
			if "Control" in file.columns[i]:
				ctrlIntensities.append(file.columns[i])
			else:
				colIntensities.append(file.columns[i])
				numIntensities += 1
	for i in range(numIntensities):
		if 'Ub' in colIntensities[i] \
		and "UbP" not in colIntensities[i]:
			file[colIntensities[i]] /= file[ctrlIntensities[0]]
		if 'UbP' in colIntensities[i]:
			file[colIntensities[i]] /= file[ctrlIntensities[1]]
		if "WCL" in colIntensities[i] \
		and "WCLP" not in colIntensities[i]:
			file[colIntensities[i]] /= file[ctrlIntensities[2]]
		if "WCLP"in colIntensities[i]:
			file[colIntensities[i]] /= file[ctrlIntensities[3]]
	whangee = []
	for i in range(numIntensities):
		if 'Whangee' in colIntensities[i]:
			whangee.append(colIntensities[i])
	Xuniques, X = np.unique(file['Protein IDs'], return_inverse=True)
	f, ax = plt.subplots(4,4)
	ax[0, 0].plot(X, file[whangee[0]], 'b.')
	ax[0, 0].set_title('Normalized' + whangee[0])
	ax[0, 1].plot(X, file_unNorm[whangee[0]], 'b.')
	ax[0, 1].set_title(whangee[0])
	
	ax[1, 0].plot(X, file[whangee[1]], 'c.')
	ax[1, 0].set_title('Normalized' + whangee[1])
	ax[1, 1].plot(X, file_unNorm[whangee[1]], 'c.')
	ax[1, 1].set_title(whangee[1])
	
	ax[2, 0].plot(X, file[whangee[2]], 'r.')
	ax[2, 0].set_title('Normalized' + whangee[2])
	ax[2, 1].plot(X, file_unNorm[whangee[2]], 'r.')
	ax[2, 1].set_title(whangee[2])
	
	ax[3, 0].plot(X, file[whangee[3]], 'm.')
	ax[3, 0].set_title('Normalized' + whangee[3])
	ax[3, 1].plot(X, file_unNorm[whangee[3]], 'm.')
	ax[3, 1].set_title(whangee[3])
	
	ax[0, 2].plot(X, file[whangee[4]], 'b.')
	ax[0, 2].set_title('Normalized' + whangee[4])
	ax[0, 3].plot(X, file_unNorm[whangee[4]], 'b.')
	ax[0, 3].set_title(whangee[4])
	
	ax[1, 2].plot(X, file[whangee[5]], 'c.')
	ax[1, 2].set_title('Normalized' + whangee[5])
	ax[1, 3].plot(X, file_unNorm[whangee[5]], 'c.')
	ax[1, 3].set_title(whangee[5])
	
	ax[2, 2].plot(X, file[whangee[6]], 'r.')
	ax[2, 2].set_title('Normalized' + whangee[6])
	ax[2, 3].plot(X, file_unNorm[whangee[6]], 'r.')
	ax[2, 3].set_title(whangee[6])
	
	ax[3, 2].plot(X, file[whangee[7]], 'm.')
	ax[3, 2].set_title('Normalized' + whangee[7])
	ax[3, 3].plot(X, file_unNorm[whangee[7]], 'm.')
	ax[3, 3].set_title(whangee[7])
	
	plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	plt.show()
	return file

def makeScatterPlot(proteinGroupsDataFrame):
    #makes scatter plot of Control vs TPK1 ko... 
    intensityList = ['Intensity Control_Ub', 'Intensity Control_UbP', 'Intensity Control_WCL', 'Intensity Control_WCLP', 'Intensity Whangee_Tpk1KO_Ub', 'Intensity Whangee_Tpk1KO_UbP', 'Intensity Whangee_Tpk1KO_WCL', 'Intensity Whangee_Tpk1KO_WCLP', 'Intensity Whangee_tunicamycin_Ub', 'Intensity Whangee_tunicamycin_UbP', 'Intensity Whangee_tunicamycin_WCL', 'Intensity Whangee_tunicamycin_WCLP']
    intenseDict = {}
    for exInt in intensityList:
        listIntensity = proteinGroupsDataFrame[exInt].tolist()
        intenseDict[exInt] = listIntensity
    f, axarr = plt.subplots(4, 2)
    axarr[0, 0].plot(intenseDict['Intensity Control_WCL'], intenseDict['Intensity Whangee_Tpk1KO_WCL'], 'b.')
    axarr[0, 1].plot(intenseDict['Intensity Control_WCLP'], intenseDict['Intensity Whangee_Tpk1KO_WCLP'], 'c.')
    axarr[1, 0].plot(intenseDict['Intensity Control_Ub'], intenseDict['Intensity Whangee_Tpk1KO_Ub'], 'r.')
    axarr[1, 1].plot(intenseDict['Intensity Control_UbP'], intenseDict['Intensity Whangee_Tpk1KO_UbP'], 'm.')

    axarr[2, 0].plot(intenseDict['Intensity Control_WCL'], intenseDict['Intensity Whangee_tunicamycin_WCL'], 'b.')
    axarr[2, 1].plot(intenseDict['Intensity Control_WCLP'], intenseDict['Intensity Whangee_tunicamycin_WCLP'], 'c.')
    axarr[3, 0].plot(intenseDict['Intensity Control_Ub'], intenseDict['Intensity Whangee_tunicamycin_Ub'], 'r.')
    axarr[3, 1].plot(intenseDict['Intensity Control_UbP'], intenseDict['Intensity Whangee_tunicamycin_UbP'], 'm.')
    plt.show()
  
  
#makeScatterPlot(proteinGroupsDataFrame)

intensity()
