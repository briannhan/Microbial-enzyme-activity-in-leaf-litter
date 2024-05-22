# -*- coding: utf-8 -*-
"""
Created on Tue May 21 22:36:11 2024

@author: Brian Chung

This script processes the Vmax data for Ashish. He intends for it to be in
the same format as the .csv file called "enzyme_all.csv". The Vmax values will
be log10 transformed, and the carbon-degrading enzymes will have their Vmax
summed up before being log10 transformed. He will then make plots with this
transformed data
"""
import os
from pathlib import Path
import pandas as pd
import numpy as np

"""Defining paths to the relevant directories and the Vmax data. Dropping the
PPO enzyme because, for some reason, he's not considering it."""
repository = Path(os.getcwd())
enzymeDataFolder = repository/"Enzyme activity data"
enzymePath = enzymeDataFolder/"Vmax.xlsx"

# Pivoting the data, log-transforming the Vmax values
indexColumns = ["timePoint", "ID", "Vegetation", "Precip"]
enzymeData = (pd.read_excel(enzymePath)
              .query("Enzyme != 'PPO'")
              .drop(columns="Replicate")
              .pivot(indexColumns, "Enzyme", "Vmax")
              .reset_index()
              )
pivotedColumns = enzymeData.columns.tolist()
colsLowercase = {column: column.lower() for column in pivotedColumns}
enzymeData = enzymeData.rename(columns=colsLowercase)
enzymeData["enzyme_c"] = enzymeData["ag"] + enzymeData["bg"] + enzymeData["bx"] + enzymeData["cbh"]
enzymes = enzymeData.columns.tolist()[4:]
for enzyme in enzymes:
    enzymeData[enzyme] = np.log10(enzymeData[enzyme])

# Converting timepoints
timePointConversion = {0: "T1", 3: "T2", 5: "T3", 6: "T4"}
for ogTimePoint in timePointConversion:
    newTimePoint = timePointConversion[ogTimePoint]
    enzymeData.loc[enzymeData.timepoint == ogTimePoint, "timepoint"] = newTimePoint

# Renaming the vegetation treatments
vegetationConversion = {"Grassland": "Grass", "CSS": "Shrub"}
for ogVeg in vegetationConversion:
    newVeg = vegetationConversion[ogVeg]
    enzymeData.loc[enzymeData.vegetation == ogVeg, "vegetation"] = newVeg

# Renaming columns
newNamesDict = {"timepoint": "time", "vegetation": "veg", "precip": "prec",
                "id": "label"}
enzymeData = enzymeData.rename(columns=newNamesDict)

# Exporting the wrangled Vmax dataset for Ashish
exportPath = enzymeDataFolder/"Vmax - log10 transformed, formatted for Ashish.csv"
if os.path.exists(exportPath) is False:
    enzymeData.to_csv(exportPath, index=False)
