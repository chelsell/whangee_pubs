import pandas
import cPickle as pic
import re
import numpy as np
import matplotlib.pyplot as plt
import urllib2
import collections
import scipy.stats as st
import os


evidenceDataFrame = pandas.DataFrame.from_csv('evidence.txt', sep='\t')
proteinGroupsDataFrame = pandas.DataFrame.from_csv('proteinGroups.txt', sep='\t')
PhosphoSTYsitesDataFrame = pandas.DataFrame.from_csv('Phospho(STY)Sites.txt', sep='\t')
peptideDataFrame = pandas.DataFrame.from_csv('peptides.txt', sep='\t')
sgdid_to_go = pic.load(open("SGDID_to_go.pkl", "rb"))
go_to_SGDID = pic.load(open("go_to_SGDID.pkl", "rb"))

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
	file = proteinGroupsDataFrame
	file_unNorm = proteinGroupsDataFrame
	sum_col = file.sum(0)
	numIntensities = 0
	colIntensities = []
	ctrlIntensities = []
	for i in range(len(file.columns)):
		if "Intensity" in file.columns[i]:
			file[file.columns[i]] /= sum_col[file.columns[i]]
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
	Xuniques, X = np.unique(file['Peptide IDs'], return_inverse=True)
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
	
def logFoldHis():
	file = intensity()
	whangee = []
	for i in range(len(file.columns)):
		if 'Whangee' in file.columns[i] \
		and 'Intensity' in file.columns[i]:
			file[file.columns[i]] = np.log2(file[file.columns[i]])
			whangee.append(file.columns[i])

	df = pandas.DataFrame({'a': file[whangee[0]], 'b': file[whangee[1]], 'c': file[whangee[2]], \
	'd': file[whangee[3]], 'e': file[whangee[4]], 'f': file[whangee[5]], \
	'g': file[whangee[6]], 'h': file[whangee[7]]})
	
	pandas.set_option('mode.use_inf_as_null', True)
	
	f, ax = plt.subplots(2,4)
	
	ax[0,0].hist(list(df.a.dropna()), bins = 30, color = 'b')
	ax[0, 0].set_title(whangee[0], fontsize = 12)
	
	ax[0, 1].hist(list(df.b.dropna()), bins = 30, color = 'c')
	ax[0, 1].set_title(whangee[1], fontsize = 12)
	
	ax[0, 2].hist(list(df.c.dropna()), bins = 30, color = 'r')
	ax[0, 2].set_title(whangee[2], fontsize = 12)
	
	ax[0, 3].hist(list(df.d.dropna()), bins = 30, color = 'm')
	ax[0, 3].set_title(whangee[3], fontsize = 12)
	
	ax[1, 0].hist(list(df.e.dropna()), bins = 30, color = 'b')
	ax[1, 0].set_title(whangee[4], fontsize = 12)
	
	ax[1, 1].hist(list(df.f.dropna()), bins = 30, color = 'c')
	ax[1, 1].set_title(whangee[5], fontsize = 12)
	
	ax[1, 2].hist(list(df.g.dropna()), bins = 30, color = 'r')
	ax[1, 2].set_title(whangee[6], fontsize = 12)
	
	ax[1, 3].hist(list(df.h.dropna()), bins = 30, color = 'm')
	ax[1, 3].set_title(whangee[7], fontsize = 12)
	plt.show()

def phosphoIntensity():
	file = PhosphoSTYsitesDataFrame
	sum_col = file.sum(0)
	numIntensities = 0
	colIntensities = []
	ctrlIntensities = []
	for i in range(len(file.columns)):
		if "Intensity" in file.columns[i] \
		and "__" not in file.columns[i]:
			file[file.columns[i]] /= sum_col[file.columns[i]]
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
	Xuniques, X = np.unique(file['id'], return_inverse=True)
	f, ax = plt.subplots(2,4)
	ax[0, 0].plot(X, file[whangee[0]], 'b.')
	ax[0, 0].set_title('Normalized' + whangee[0])
	
	ax[0, 1].plot(X, file[whangee[1]], 'c.')
	ax[0, 1].set_title('Normalized' + whangee[1])
	
	ax[0, 2].plot(X, file[whangee[2]], 'r.')
	ax[0, 2].set_title('Normalized' + whangee[2])
	
	ax[0, 3].plot(X, file[whangee[3]], 'm.')
	ax[0, 3].set_title('Normalized' + whangee[3])
	
	ax[1, 0].plot(X, file[whangee[4]], 'b.')
	ax[1, 0].set_title('Normalized' + whangee[4])
	
	ax[1, 1].plot(X, file[whangee[5]], 'c.')
	ax[1, 1].set_title('Normalized' + whangee[5])
	
	ax[1, 2].plot(X, file[whangee[6]], 'r.')
	ax[1, 2].set_title('Normalized' + whangee[6])
	
	ax[1, 3].plot(X, file[whangee[7]], 'm.')
	ax[1, 3].set_title('Normalized' + whangee[7])
	return file

def phosphoLogFoldHis():
	file = phosphoIntensity()
	whangee = []
	for i in range(len(file.columns)):
		if 'Whangee' in file.columns[i]:
			if 'Intensity' in file.columns[i] and '__' not in file.columns[i]:
				file[file.columns[i]] = np.log2(file[file.columns[i]])
				whangee.append(file.columns[i])
	df = pandas.DataFrame({'a': file[whangee[0]], 'b': file[whangee[1]], 'c': file[whangee[2]], \
	'd': file[whangee[3]], 'e': file[whangee[4]], 'f': file[whangee[5]], \
	'g': file[whangee[6]], 'h': file[whangee[7]]})
	
	pandas.set_option('mode.use_inf_as_null', True)
	f, ax = plt.subplots(2,4)

	ax[0,0].hist(list(df.a.dropna()), color = 'b')
	ax[0, 0].set_title(whangee[0], fontsize = 12)
	
	ax[0, 1].hist(list(df.b.dropna()), color = 'c')
	ax[0, 1].set_title(whangee[1], fontsize = 12)
	
	ax[0, 2].hist(list(df.c.dropna()), color = 'r')
	ax[0, 2].set_title(whangee[2], fontsize = 12)
	
	ax[0, 3].hist(list(df.d.dropna()), color = 'm')
	ax[0, 3].set_title(whangee[3], fontsize = 12)
	
	ax[1, 0].hist(list(df.e.dropna()), color = 'b')
	ax[1, 0].set_title(whangee[4], fontsize = 12)
	
	ax[1, 1].hist(list(df.f.dropna()), color = 'c')
	ax[1, 1].set_title(whangee[5], fontsize = 12)
	
	ax[1, 2].hist(list(df.g.dropna()), color = 'r')
	ax[1, 2].set_title(whangee[6], fontsize = 12)
	
	ax[1, 3].hist(list(df.h.dropna()), color = 'm')
	ax[1, 3].set_title(whangee[7], fontsize = 12)
	plt.show()

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

def calculate_enrichment(intensityStringName, proteinGroupsDataFrame, go_to_SGDID,sgdidList, n=100):

    intensityList = list(proteinGroupsDataFrame[intensityStringName])
    gene_data = []
    for index in xrange(0,len(intensityList)):
        gene_data.append((sgdidList[index],intensityList[index]))

    # first, sort the gene data by expression value, from top to bottom
    gene_list = sorted(gene_data, key=lambda g: g[1], reverse=True)

    # gene_set1 is a set of the N highest-expressing genes
    gene_set1 = set(g for g,v in gene_list[:n])

    # gene_set2 is a set of the N lowest-expressing genes
    gene_set2 = set(g for g,v in gene_list[-n:])

    # dictionaries that map from goid -> -log p-value
    e_scores = dict() # positive-enrichment scores
    ne_scores = dict() # negative-enrichment scores

    # for each goid in the go_to_genes dictionary...
    for go in go_to_SGDID:
        # convert the list of genes into a set
        goset = set(go_to_SGDID[go])

        # The function's interface is st.hypergeom.logsf(x, M, n, N),
        # and it returns ln(p-value) so I convert to log10.

        # I'm using set intersection to figure out how big x is.
        # I substract one because the sf calculates "more than x"
        # and I want "x or more" (e.g. "more than x-1")
        e_scores[go] = -(st.hypergeom.logsf(len(goset & gene_set1) - 1,
                                            len(gene_list),
                                            len(goset),
                                            len(gene_set1)) / np.log(10))

        # I do the same thing for the negative enrichment, just using the
        # genes from the bottom of the list this time
        ne_scores[go] = -(st.hypergeom.logsf(len(goset & gene_set2) - 1,
                                             len(gene_list),
                                             len(goset),
                                             len(gene_set2)) / np.log(10))


    # I convert the dictionaries into lists by sorting them by their scores.
    e_scores = sorted(((go,e_scores[go]) for go in e_scores),
                      key=lambda (go,s): s, reverse=True)
    ne_scores = sorted(((go,ne_scores[go]) for go in ne_scores),
                       key=lambda (go,s): s, reverse=True)
    #print e_scores, ne_scores
    return e_scores,ne_scores


def checkGoLengths(control, test):
    #print len(e0), len(e1)
    gocontrol = []
    for e in control:
        gocontrol.append(e[0])
    gotest = []
    for e in test:
        gotest.append(e[0])
    controlnottest = list(set(gocontrol)-set(gotest))
    for e in controlnottest:
        test.append((e, -0.0))
    testnotcontrol = list(set(gotest)-set(gocontrol))
    for e in testnotcontrol:
        control.append((e,-0.0))
    maxlist = [control[0][1], test[0][1]]
    test_scores = []
    goIDlist = []
    for t in sorted(test, key=lambda g: g[0], reverse=True):
        test_scores.append(t[1])
        goIDlist.append(t[0])
    control_scores = []
    for c in sorted(control, key=lambda g: g[0], reverse=True):
        control_scores.append(c[1])
    return test_scores, control_scores, goIDlist, maxlist



def printEnrichment(controlVsTestTupleList, proteinGroupsDataFrame, go_to_SGDID,sgdidList):
    f, axarr = plt.subplots(4, 2)
    f.subplots_adjust(hspace = 0.75)
    i = 0
    j = 0
    font ={'family':'normal', 'weight':'bold', 'size':8}
    plt.rc('font', **font)

    for (control, test) in controlVsTestTupleList:
        e1, n1 = calculate_enrichment(test, proteinGroupsDataFrame, go_to_SGDID,sgdidList)
        e0, n0 = calculate_enrichment(control, proteinGroupsDataFrame, go_to_SGDID,sgdidList)
        e0list, e1list, egoidlist, emaxlist = checkGoLengths(e0, e1)
        n0list, n1list, ngoidlist, nmaxlist = checkGoLengths(n0, n1)
        maxlist = [emaxlist[0], emaxlist[1], nmaxlist[0], nmaxlist[1]]
        axarr[i, j].plot(e0list, e1list, 'b.')
        axarr[i, j].plot(n0list, n1list, 'r.')
        axarr[i, j].set_title(control +'vs'+ test, fontsize=8)
        axarr[i, j].set_xlabel('control')
        axarr[i, j].set_ylabel('test') 
        axarr[i, j].set_xlim([0, max(maxlist)])
        axarr[i, j].set_ylim([0, max(maxlist)])
        """index = 0
        #to add all go ids to plot
        for xy in zip(e0list, e1list):                                               
            axarr[i, j].annotate('%s' % egoidlist[index], xy=xy, textcoords='data')
        index+=1"""
        if j == 1:
            i+=1
            j=0
        else:
            j+=1
    plt.show()


def getGo(sgdid_to_go, proteinGroupsDataFrame, go_to_SGDID, controlVsTestTupleList):
    """Add go ids to DataFrame
    """
    fastaList = list(proteinGroupsDataFrame['Fasta headers'])
    sgdidList = []
    for header in fastaList:
        sgdid = ("".join(re.findall(r'SGDID:S[0-9]*', header))).split(':')
        if len(sgdid)>1:
            sgdidList.append(sgdid[1])
        else:
            sgdidList.append(0)
    golist = []
    for sgdid in sgdidList:
        if sgdid_to_go[sgdid]:
            golist.append("|".join(sgdid_to_go[sgdid]))
        else:
            golist.append(0)
    proteinGroupsDataFrame['goid'] = golist
    printEnrichment(controlVsTestTupleList, proteinGroupsDataFrame, go_to_SGDID,sgdidList)
  
#makeScatterPlot(proteinGroupsDataFrame)

intensity()
