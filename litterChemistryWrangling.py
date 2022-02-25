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
functional = functionalOG.drop(columns=["carbo1", "carbo3", "amide1",
                                        "amide2"])
# %%
# Renaming independent variable columns and the values in these columns
renameColsDict = {"veg": "Vegetation", "prec": "Precip", "time": "timePoint"}
functional.rename(renameColsDict, axis="columns", inplace=True)

# Renaming the values in each independent variable column
for index, row in functional.iterrows():
    # Renaming vegetation values
    if row.Vegetation == "Grass":
        functional.loc[index, "Vegetation"] = "Grassland"

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
# Exporting the data. I'm wondering if I have to analyze lignin later, so just
# to clarify, this data only contains spectral area of carbohydrates and
# proteins
fileName = "Carbohydrates and Proteins FTIR.xlsx"
filePath = litterChemFolderPath/fileName
functional.to_excel(filePath, index=False)
