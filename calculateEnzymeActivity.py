# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 23:24:06 2020

@author: Brian Chung
The purpose of this script is to calculate the enzyme activity of the litter
from Loma Ridge. This script processes enzyme activity & litter dry weight data
"""
# %%
# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as py
import os
from pathlib import Path

# %%
# Making file & folder paths to read in data
workingDirectory = Path(os.getcwd())
enzymeFolder = "Enzyme activity data"
dryMassFolder = "Litter dry weights"
enzymeFolderPath = workingDirectory/enzymeFolder
enzymeFiles = os.listdir(enzymeFolderPath)
dryMassFolderPath = workingDirectory/dryMassFolder
dryMassFiles = os.listdir(dryMassFolderPath)

# %%
# Processing dry weight data before using it for enzyme activity calculations
# Processing includes (1) renaming columns; (2) Filling in time points;
# (3) Adding columns for wet and dry assay mass and proportion of leaf litter
# as dry mass; (4) Calculating proportion dry assay mass; and (5) isolating
# timepoints T0, T3, and T5 (T6 will be isolated from a different spreadsheet)

# The following line makes pandas print all columns
pd.set_option("display.max_columns", None)

# Reading in dry mass file
dryMassFilePath = dryMassFolderPath/dryMassFiles[0]
dryDF = pd.read_excel(dryMassFilePath)

# (1) Renaming columns
oldColumns = dryDF.columns.tolist()
newColumns = ["Envelope mass (g)", "Envelope + wet (g)",
              "Envelope + dry (g)", "Wet litter mass (g)",
              "Dry litter mass (g)"]
newColumnsDict = {oldColumns[2]: newColumns[0], oldColumns[3]: newColumns[1],
                  oldColumns[4]: newColumns[2], oldColumns[5]: newColumns[3],
                  oldColumns[6]: newColumns[4]}
dryDF = dryDF.rename(columns=newColumnsDict)


# (2) Filling the time column with time points
for index, row in dryDF.iterrows():
    currentTimepointType = type(row["Time"])
    if currentTimepointType != str:
        previousRow = dryDF.iloc[index - 1]
        dryDF.loc[index, "Time"] = previousRow["Time"]

# (3) & (4) Adding columns & calculating dry mass
dryDF["Wet assay (g)"] = 0.4
dryDF["Dry proportion (%)"] = 100*(dryDF["Dry litter mass (g)"] /
                                   dryDF["Wet litter mass (g)"])
dryDF["Dry assay (g)"] = dryDF["Wet assay (g)"]*dryDF["Dry proportion (%)"]/100

# (5) Isolating timepoints T0, T3, T5
timepoints = dryDF.groupby("Time")["Time"].count().index.tolist()
timepoints = [timepoints[0], timepoints[3], timepoints[5]]
T0_3_5Bool = dryDF["Time"].isin(timepoints)
dryDF = dryDF[T0_3_5Bool]

# I'm going to test with a single enzyme file & dry weights file
# Reading in 1 enzyme file and the dry weights excel that contains dry weights
# for all time points
