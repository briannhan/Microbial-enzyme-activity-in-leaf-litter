# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:29:55 2022

@author: Brian Chung
This script tests whether the Box-Cox transformation can be used in multiple
comparisons, and will use the transformation from scipy.stats. As the
transformation transforms a whole dataset across all treatments, I want to see
if this transformation will remove the patterns that have been observed in the
untransformed data. This concerns stems from me thinking that this
transformation, when done, will transform different treatments so that they
are from the same population, which would mean that parametric inference would
not detect a significant difference between treatments.

I will test the transformation on 2 functional groups in the litter chemistry
data: the FTIR spectral area for (1) amide1 and (2) glycosidic bonds. I will
perform Tukey's Honest Significant Difference tests afterwards. The patterns
I seek are:

- a vegetation x precipitation interaction for glycosidic bonds in which
grassland x ambient > grassland x drought > CSS ambient = CSS drought
- a timepoint main effect where the amide1 increases with time
- a precipitation effect where amide1 is higher in drought
"""

import os
from pathlib import Path
import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import MultiComparison as MC
import cld

# Reading in the data
cwd = Path(os.getcwd())
cwdFilesNfolders = os.listdir(cwd)

cwdFolders = []
for item in cwdFilesNfolders:
    itemPath = cwd/item
    if os.path.isdir(itemPath):
        cwdFolders.append(item)

litterChemFolder = cwd/'Litter chemistry'
litterChemFiles = os.listdir(litterChemFolder)
ogDataPath = litterChemFolder/'Carbohydrates and Proteins FTIR.xlsx'
ogData = pd.read_excel(ogDataPath)
ogData = ogData[(ogData.functionalGroup == "amide1")
                | (ogData.functionalGroup == "glycosidicBond")]
ogData["timePoint"] = ogData["timePoint"].astype(str)

# Transforming the data
dataTransformed = ogData.copy()
amide1 = dataTransformed.loc[dataTransformed.functionalGroup == "amide1", "value"]
amide1trans, amide1lambda = stats.boxcox(amide1)
dataTransformed.loc[dataTransformed.functionalGroup == "amide1", "value"] = amide1trans
glycosidic = dataTransformed.loc[dataTransformed.functionalGroup == "glycosidicBond", "value"]
glycosidicTrans, glycosidicLambda = stats.boxcox(glycosidic)
dataTransformed.loc[dataTransformed.functionalGroup == "glycosidicBond", "value"] = glycosidicTrans
# %%
# Purpose: Performing Tukey's HSD test

# (1) Obtaining the folder path that will contain the raw Tukey results for
# this script
statAnalysesFolder = cwd/'Statistical analyses'
statAnalysesFilesNfolders = os.listdir(statAnalysesFolder)
normalityFolder = statAnalysesFolder/'Normality testing'
normalityFilesNfolders = os.listdir(normalityFolder)
boxcoxFolder = normalityFolder/'Box Cox testing'


# (2) Writing a function to translate Tukey's results from a statsmodel
# SimpleTable object over to a pandas dataframe
def translateResults(SimpleTable, fileName):
    """
    Translates a statsmodels SimpleTable object that contains Tukey results,
    over to a Pandas dataframe. This will be done by manually reading and
    writing each line in the SimpleTable over to a text file, which will be
    read in using the read_csv() function from pandas to create a dataframe

    Parameters
    ----------
    SimpleTable : statsmodels SimpleTable
        This object is native to statsmodels and contains the raw Tukey
        results that will be translated over to a Pandas dataframe.
    fileName : str
        The name of the text file that will be created. Does not end with
        .txt

    Returns
    -------
    Pandas dataframe that contains the raw Tukey results, having been
    translated from a statsmodels SimpleTable format.

    """
    fileName = fileName + ".txt"
    filePath = boxcoxFolder/fileName
    file = open(filePath, "w")
    for row in SimpleTable:
        for element in row:
            element = str(element)
            file.write(element)
            file.write(',')
        file.write('\n')
    file.close()
    resultsDF = pd.read_csv(filePath, header=0)
    return resultsDF


# (3) Looking at Tukey's results of time on transformed amide1
amide1DF = dataTransformed[dataTransformed.functionalGroup == "amide1"]
amide1timeMC = MC(amide1DF.value, amide1DF.timePoint)
amide1timeST = amide1timeMC.tukeyhsd().summary()  # ST = SimpleTable format
amide1timeDF = translateResults(amide1timeST, "amide1, time")
amide1timeCLD = cld.main(amide1timeDF)
"""Ok results match Tukey on untransformed amide1"""

# (4) Performing Tukey's HSD on precipitation on transformed amide1
amide1pptMC = MC(amide1DF.value, amide1DF.Precip)
amide1pptST = amide1pptMC.tukeyhsd().summary()
amide1pptDF = translateResults(amide1pptST, "amide1, precipitation")
amide1pptCLD = cld.main(amide1pptDF)
"""Ok this matches Tukey's of precipitation on untransformed amide1"""

# (5) Performing Tukey's HSD on vegetation x precipitation on transformed
# glycosidic bond FTIR spectral area
glycosidicDF = dataTransformed[dataTransformed.functionalGroup == "glycosidicBond"]
timeXvegSeries = glycosidicDF.Vegetation + " x " + glycosidicDF.Precip
glycosidicMC = MC(glycosidicDF.value, timeXvegSeries)
glycosidicST = glycosidicMC.tukeyhsd().summary()
glycosidicDF = translateResults(glycosidicST, "glycosidic, veg x ppt")
glycosidicCLD = cld.main(glycosidicDF)
"""And this also matches Tukey's of vegetation x precipitation on untransformed
glycosidic bonds. So, in conclusion, the Box-Cox transformation seems like a
valid transformation to do for multiple comparisons. It will not transform
different treatment groups into a single population."""
