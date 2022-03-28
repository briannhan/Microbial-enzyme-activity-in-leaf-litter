# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 16:43:16 2022

@author: Brian Chung
This script wrangles the CAZymes gene domain relative abundance to fit the
existing code I use to check data normality and perform ANOVAs and Tukey's.
The original dataset, as Ashish had given me, consists of columns that describe
- sample ID
- vegetation type (independent variable, blocking factor)
- precipitation (independent variable)
- time point (independent variable)
- the relative abundance of gene domains that produce CAZymes that degrade a
particular substrate
- a column describing the total number of gene domains that produce CAZymes

The relative abundance for a particular substrate, in this case, is defined as

number of gene domains for a particular substrate
-------------------------------------------------
total number of CAZyme gene domains
"""
import pandas as pd
import os
from pathlib import Path

cwd = Path(os.getcwd())
cwdFilesNfolders = os.listdir(cwd)

# Obtaining the names of folders in the working directory
cwdFolders = []
for item in cwdFilesNfolders:
    itemPath = cwd/item
    if os.path.isdir(itemPath):
        cwdFolders.append(item)

# Reading in the CAZymes gene domains data
CAZymeFolder = cwd/'CAZyme metagenomic data'
CAZymeFolderContents = os.listdir(CAZymeFolder)
CAZymePath = CAZymeFolder/'cazymes by substrate_final.csv'
CAZymeData = pd.read_csv(CAZymePath)
# %%
# Purpose: Wrangling the CAZyme data

# (1) Melting the various substrate columns down into fewer columns so I can
# display all columns in the console. The "value" column will be the dependent
# variable column
idVars = ["id", "veg", "prec", "samp"]
CAZymeData = CAZymeData.melt(idVars, var_name="Substrate")

# (2) Creating a parameter column to describe whether the dependent variable is
# the relative abundance for a particular substrate or the total number of
# CAZyme gene domains
substrates = set(CAZymeData.Substrate.tolist())
substrates = list(substrates)
if "Total" in substrates:
    substrates.remove("Total")
for substrate in substrates:
    boolean = CAZymeData["Substrate"] == substrate
    CAZymeData.loc[boolean, "Parameter"] = "Relative abundance"
totalBoolean = CAZymeData["Substrate"] == "Total"
CAZymeData.loc[totalBoolean, "Parameter"] = "Total CAZyme domains"

# (3) Renaming the sample ID column, the independent variable columns, and
# their values so
# that their new names match the names I use for them in my existing normality,
# factorial ANOVA, and Tukey code
renameColumns = {"id": "ID", "samp": "timePoint", "veg": "Vegetation",
                 "prec": "Precip"}
CAZymeData.rename(columns=renameColumns, inplace=True)
timePointOGvals = set(CAZymeData.timePoint.tolist())
for timePoint in timePointOGvals:
    timePointBool = CAZymeData["timePoint"] == timePoint
    if timePoint == "T1":
        CAZymeData.loc[timePointBool, "timePoint"] = "0"
    elif timePoint == "T2":
        CAZymeData.loc[timePointBool, "timePoint"] = "3"
    elif timePoint == "T3":
        CAZymeData.loc[timePointBool, "timePoint"] = "5"
    elif timePoint == "T4":
        CAZymeData.loc[timePointBool, "timePoint"] = "6"
vegOGvals = set(CAZymeData.Vegetation.tolist())
for vegetation in vegOGvals:
    vegBool = CAZymeData["Vegetation"] == vegetation
    if vegetation == "grass ":
        CAZymeData.loc[vegBool, "Vegetation"] = "Grassland"
    elif vegetation == "shrub ":
        CAZymeData.loc[vegBool, "Vegetation"] = "CSS"
pptOGvals = set(CAZymeData.Precip.tolist())
for pptTreatment in pptOGvals:
    pptBool = CAZymeData["Precip"] == pptTreatment
    if pptTreatment == "amb ":
        CAZymeData.loc[pptBool, "Precip"] = "Ambient"
    elif pptTreatment == "red ":
        CAZymeData.loc[pptBool, "Precip"] = "Reduced"

# (4) Exporting the fully wrangled CAZyme gene domain data
exportPath = CAZymeFolder/"Wrangled CAZyme gene domains.xlsx"
if os.path.exists(exportPath) is False:
    CAZymeData.to_excel(exportPath, index=False)
