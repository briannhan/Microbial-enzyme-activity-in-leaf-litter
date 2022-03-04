# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:23:48 2022

@author: Brian Chung
This script checks the compact letter displays (CLDs) that I manually created
for enzyme Vmax and Km when I was working on my thesis to see if my work by
hand matches up with the raw Tukey test results I used to create these CLDs.
"""
import pandas as pd
import cld
import os
from pathlib import Path

cwd = Path(os.getcwd())
filesNfolders = os.listdir()
cwdFolders = []
cwdFiles = []
for name in filesNfolders:
    namePath = cwd/name
    if os.path.isdir(namePath):
        cwdFolders.append(name)
    elif os.path.isfile(namePath):
        cwdFiles.append(name)

statsFolderPath = cwd/'Statistical analyses'
tukeyFolderPath = statsFolderPath/"Tukey posthoc"
enzymes = ["BG", "BX", "CBH", "LAP", "NAG", "PPO"]
indeVar = ["timePoint", "Precipitation", "Vegetation",
           "timePoint x Precipitation", "timePoint x Vegetation",
           "Vegetation x Precipitation", "Three-way"]

# Start the process of checking compact letter displays that I manually created
# for Vmax and Km of all enzymes
for enzyme in enzymes:
    enzymeResultsFolder = tukeyFolderPath/enzyme
    resultsFiles = os.listdir(enzymeResultsFolder)

    # Looking for CLD files in each enzyme folder as well as raw Tukey results
    excelFiles = []
    for file in resultsFiles:
        if file.endswith(".xlsx") and not file.endswith("groups.xlsx"):
            excelFiles.append(file)

    # Dividing the Excel files into raw Tukey results and
    # CLDs created from the raw results
    strippedExcelFiles = []
    for file in excelFiles:
        file = file.strip(".xlsx")
        strippedExcelFiles.append(file)
    tukeyResultsFiles = []
    cldFiles = []
    for file in excelFiles:
        if "annotated" in file:
            cldFiles.append(file)
        elif "annotated" not in file:
            tukeyResultsFiles.append(file)
    listOfLists = []  # Will be populated with lists where each list contains
    # a pair of files, a Tukey results file and its associated cld file
    for i in range(len(tukeyResultsFiles)):
        for j in range(len(cldFiles)):
            tukeyFile = tukeyResultsFiles[i]
            cldFile = cldFiles[j]
            if tukeyFile in cldFile:
                innerList = [tukeyFile, cldFile]
                listOfLists.append(innerList)

    # Now, actually checking to see if cld files match up with Tukey results
    # files
    for pair in listOfLists:
        tukeyFileName = pair[0] + ".xlsx"
        tukeyFilePath = enzymeResultsFolder/tukeyFileName
        tukeyFile = pd.read_excel(tukeyFilePath)
        cldFileName = pair[1] + ".xlsx"
        cldFilePath = enzymeResultsFolder/cldFileName
        cldFile = pd.read_excel(cldFilePath)
        result = cld.checkCLD(tukeyFile, cldFile)
        if result == "Bad":
            print(enzyme)
            print(pair)
            print('\n')
"""Ok, I don't see anything printed out, so this means that my manual compact
letter displays are actually good. What do you know."""
