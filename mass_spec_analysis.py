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
  
  
makeScatterPlot(proteinGroupsDataFrame)
