# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 13:58:54 2022

@author: Brian Chung
This script wrangles FTIR litter chemistry data into a format that is more
conducive for analysis. Modifications to the data by this script includes, but
not in any particular order:

(1) Combining the 2 amide columns into 1 column and the carbohydrate esters
columns into 1 column. Columns will be combined by simply adding them together
(2) Renaming the independent variable columns and the values in these columns

Afterwards, factorial ANOVAs and Tukey's post-hoc will be conducted on this
data.
"""

import pandas as pd
from pathlib import Path
import os

# Reading in the data
cwd = Path(os.getcwd())
litterChemFolderPath = cwd/"Litter chemistry"
litterChemFiles = os.listdir(litterChemFolderPath)
litterChemPath = litterChemFolderPath/'ftir.percent.area_brian.xlsx'
functionalOG = pd.read_excel(litterChemPath)
# %%
# Combining amide and carbohydrate ester columns
functionalOG["carboEster"] = functionalOG.carbo1 + functionalOG.carbo3
functionalOG["amide"] = functionalOG.amide1 + functionalOG.amide2

# Dropping the original carbohydrate ester and amide columns
# functional = functionalOG.drop(columns=["carbo1", "carbo3", "amide1",
#                                         "amide2"])
# %%
# Renaming independent variable columns and the values in these columns
renameColsDict = {"veg": "Vegetation", "prec": "Precip", "time": "timePoint",
                  "carbo1": "carboEster1", "carbo2": "glycosidicBond",
                  "carbo3": "carboEster2", "carbo4": "C_O_stretching"}
# functional.rename(renameColsDict, axis="columns", inplace=True)
functional = functionalOG.rename(renameColsDict, axis="columns")

# Renaming the values in each independent variable column
for index, row in functional.iterrows():
    # Renaming vegetation values
    if row.Vegetation == "Grass":
        functional.loc[index, "Vegetation"] = "Grassland"
    elif row.Vegetation == "Shrub":
        functional.loc[index, "Vegetation"] = "CSS"

    # Renaming time points
    if row.timePoint == "T1":
        functional.loc[index, "timePoint"] = "0"
    elif row.timePoint == "T2":
        functional.loc[index, "timePoint"] = "3"
    elif row.timePoint == "T3":
        functional.loc[index, "timePoint"] = "5"
    elif row.timePoint == "T4":
        functional.loc[index, "timePoint"] = "6"
# %%
# Purpose: Calculate FTIR spectral area of different raw FTIR files to see
# which files are actually the raw FTIR data that Ashish used to calculate
# the carbo1, carbo2, carbo3, carbo4, amide1, amide2, lipid, and alkane
# spectral areas that he then sent me

# (1) Calculating carbo1 spectral area for sample at plot 14 and timepoint 0
# using the ftir.txt file. This spectral area has a range of 970-1015 cm-1
ftirPath = litterChemFolderPath/"ftir.txt"
ftirTxt = pd.read_csv(ftirPath, sep="\t", header=0)
ftirTxt = ftirTxt[(ftirTxt.wavelength >= 970) & (ftirTxt.wavelength <= 1015)]
ftirTxt14 = ftirTxt["t0-14"].sum()
print("ftirTxt14:", ftirTxt14)
print('\n')

# (2) Doing the same but instead of calculating using ftir.txt, I calculate
# using ftir.norm.csv
ftirNormPath = litterChemFolderPath/"ftir.norm.csv"
ftirNorm = pd.read_csv(ftirNormPath)
ftirNormOGcols = ftirNorm.columns.tolist()
ftirNormWavelengths = ftirNormOGcols[1:]
ftirNorm = ftirNorm.melt(id_vars="id", value_vars=ftirNormWavelengths,
                         var_name="wavelength")
ftirNorm["wavelength"] = ftirNorm["wavelength"].str.lstrip("x")
ftirNorm["wavelength"] = ftirNorm["wavelength"].astype(float)
ftirNorm = ftirNorm[(ftirNorm["wavelength"] >= 970)
                    & (ftirNorm["wavelength"] <= 1015)
                    & (ftirNorm.id == "t0-14")]
ftirNorm14 = ftirNorm["value"].sum()
print("ftirNorm14:", ftirNorm14)
print('\n')
"""So none of these 2 numbers match the value that Ashish gave me for carbo1
of plot plot 14 at time T0, which is about 7.768% of FTIR spectral area. So, I
don't know how he calculated these numbers that he gave me in the first place.
"""
# %%
# Exporting the data. I'm wondering if I have to analyze lignin later, so just
# to clarify, this data only contains spectral area of carbohydrates and
# proteins

# Collapsing the functional group columns into 2 columns: the functional group
# and the spectral area occupied by the functional group
functionalGroupsOG = ["glycosidicBond", "C_O_stretching", "alkane", "lipid",
                      "carboEster", "amide", "carboEster1", "carboEster2",
                      "amide1", "amide2"]
idVars = ["id", "Vegetation", "Precip", "timePoint"]
functional = functional.melt(idVars, functionalGroupsOG, "functionalGroup")

# Exporting the newly wrangled dataframe
fileName = "Carbohydrates and Proteins FTIR.xlsx"
filePath = litterChemFolderPath/fileName
functional.to_excel(filePath, index=False)
