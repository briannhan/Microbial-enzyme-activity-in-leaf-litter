# -*- coding: utf-8 -*-
"""
Created on Sat May 29 14:55:28 2021

@author: Brian Chung
This script performs linear regressions on the Vmax and Km of each enzyme to
test hypothesis 1, which states that amount of enzymes determine the amount of
Km in such a way that both are positively correlated.
"""

import pandas as pd
import os
from pathlib import Path
from scipy.stats import linregress
import matplotlib.pyplot as py
from matplotlib.patches import Patch
import numpy as np
from datetime import datetime
start = datetime.now()

# Reading in parameters file
cwd = Path(os.getcwd())
cwdDirsFiles = os.listdir(cwd)
activityFolder = cwd/'Enzyme activity data'
activityContents = os.listdir(activityFolder)
paramPath = activityFolder/'Parameters - log 10 transformed.xlsx'
parameters = pd.read_excel(paramPath)
paramCols = parameters.columns.tolist()
statAnalysesFolder = cwd/'Statistical analyses'
statAnalysesContents = os.listdir(statAnalysesFolder)
linRegressFolder = statAnalysesFolder/'Vmax-Km linear regressions'
enzymes = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG", "PPO"]
pList = []
rList = []
r2List = []
# %%
# Purpose: Performing linear regression on AG as a proof of concept

# (1) Performing linear regression on AG
AGparams = parameters[parameters.Enzyme == "AG"]
VmaxAG = AGparams[AGparams.Parameter == "Vmax"]["value"].tolist()
KmAG = AGparams[AGparams.Parameter == "Km"]["value"].tolist()
resultsAG = linregress(VmaxAG, KmAG)
rAG = resultsAG.rvalue
rAGstr = "R = {0:.3f}".format(rAG)
rList.append(rAG)
r2AG = rAG**2
r2AGstr = "R-squared = {0:.3f}".format(r2AG)
r2List.append(r2AG)
pAG = resultsAG.pvalue
pAGstr = "p = {0:.3f}".format(pAG)
pList.append(pAG)

# (2) Plotting Km vs Vmax for AG
slopeAG = resultsAG.slope
bAG = resultsAG.intercept
maxVmaxAG = max(VmaxAG)
minVmaxAG = min(VmaxAG)
generatedVmaxAG = np.linspace(minVmaxAG, maxVmaxAG)
KmAGestimated = (slopeAG*generatedVmaxAG) + bAG
colorEstimatedVals = (170/255, 51/255, 119/255)
colorActualVals = (68/255, 119/255, 170/255)
py.figure("Vmax-Km linear regressions", (20, 14))
AGplotPos = enzymes.index("AG") + 1
py.subplot(2, 4, AGplotPos)
py.title("AG")
yAxis = r"$log_{10}$ $(K_{m})$"
py.ylabel(yAxis)
xAxis = r"$log_{10}$ $(V_{max})$"
py.xlabel(xAxis)
actualDataAG = py.plot(VmaxAG, KmAG, marker=".",
                       markerfacecolor=colorActualVals, linestyle="None")
AGregress = py.plot(generatedVmaxAG, KmAGestimated,
                    markerfacecolor=colorEstimatedVals, linestyle="-")
rPatch = Patch(color="w")
r2Patch = Patch(color="w")
pPatch = Patch(color="w")
legendHandles = [actualDataAG[0], AGregress[0], pPatch, r2Patch, rPatch]
legendLabels = ["Actual data", "Regression", pAGstr, r2AGstr, rAGstr]
py.legend(legendHandles, legendLabels)
# %%
# Purpose: Performing linear regression for all enzymes

# (1) Clearing everything to make a new plot
pList = []
rList = []
r2List = []
py.close()
py.figure("Vmax-Km linear regressions", (20, 14))


# (2) Abstracting the linear regressions and plotting process
def linearRegression():
    """
    Performs Vmax-Km linear regressions for every enzyme and plots the actual
    data and trend line for each enzyme. Also updates the lists of p-values,
    R, and R-squared values.

    Returns
    -------
    None.

    """
    for enzyme in enzymes:
        ezParams = parameters[parameters.Enzyme == enzyme]
        Vmax = ezParams[ezParams.Parameter == "Vmax"]["value"].tolist()
        Km = ezParams[ezParams.Parameter == "Km"]["value"].tolist()
        results = linregress(Vmax, Km)
        r = results.rvalue
        rStr = "R = {0:.3f}".format(r)
        rList.append(r)
        r2 = r**2
        r2Str = "R-squared = {0:.3f}".format(r2)
        r2List.append(r2)
        p = results.pvalue
        pList.append(p)
        if p < 0.05 and p >= 0.01:
            pStr = "p < 0.05"
        elif p < 0.01 and p >= 0.001:
            pStr = "p < 0.01"
        elif p < 0.001:
            pStr = "p < 0.001"

        # Plotting Km vs Vmax for AG
        m = results.slope
        b = results.intercept
        maxVmax = max(Vmax)
        minVmax = min(Vmax)
        generatedVmax = np.linspace(minVmax, maxVmax)
        KmEstimated = (m*generatedVmax) + b
        subplotPos = enzymes.index(enzyme) + 1
        py.subplot(2, 4, subplotPos)
        py.title(enzyme)
        yAxis = r"$log_{10}$ $(K_{m})$"
        py.ylabel(yAxis)
        xAxis = r"$log_{10}$ $(V_{max})$"
        py.xlabel(xAxis)
        actualData = py.plot(Vmax, Km, marker=".",
                             markerfacecolor=colorActualVals, linestyle="None")
        regress = py.plot(generatedVmax, KmEstimated,
                          markerfacecolor=colorEstimatedVals, linestyle="-")
        legendHandles = [actualData[0], regress[0], pPatch, r2Patch, rPatch]
        legendLabels = ["Actual data", "Regression", pStr, r2Str, rStr]
        py.legend(legendHandles, legendLabels)
    figPath = linRegressFolder/"Vmax-Km linear regressions.png"
    if not os.path.exists(figPath):
        py.savefig(figPath)
    return


# (3) Performing & plotting linear regressions
linearRegression()

# (5) Saving p-values, R, and R-squared values as an Excel file
resultsDF = pd.DataFrame({"Enzyme": enzymes, "pValue": pList, "R": rList,
                          "R2": r2List})
exportName = "Linear regression results.xlsx"
exportPath = linRegressFolder/exportName
resultsDF.to_excel(exportPath, index=False)
# %%
print(datetime.now() - start)
