# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 18:51:39 2020

@author: Brian Nhan Thien Chung
This is a module that is meant to facilitate the wrangling of enzyme activity
data and to help plot them.
"""
import pandas as pd

# %%
# Functions to wrangle both hydrolytic & oxidative enzyme data


def processLongNamesAndWells(enzymeData):
    """Processes the long sample names & wells in an enzyme dataframe that
    occurred from reading in a dataframe of raw enzyme data. The long sample
    names are in a column, which are then split up into 2 columns. The Well
    column will also be used to create 2 new columns


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of raw enzyme data from reading in a raw enzyme data text
        file. Has 3 columns: Well, Long sample name, Assay. Well will be used
        to create 2 new columns, while Long sample name will also be used
        to create 2 new columns

    Returns
    -------
    The dataframe with the names & wells processed & with new columns. Now,
    the dataframe can be checked for missing, extra, misnamed plates

    """
    enzymeData["Long sample name"] = enzymeData["Long sample name"].str.split("_")
    for index, row in enzymeData.iterrows():
        longSampleName = row["Long sample name"]
        enzymeData.loc[index, "Assay date"] = longSampleName[0]
        if "B" in longSampleName:
            enzymeData.loc[index, "ID"] = "B"
        elif "X" in longSampleName[2]:
            enzymeData.loc[index, "ID"] = longSampleName[2]
        well = row["Well"]
        wellRow = well[0]
        wellColumn = int(well[1:])
        enzymeData.loc[index, "PlateRow"] = wellRow
        enzymeData.loc[index, "PlateCol"] = wellColumn
    enzymeData = enzymeData.drop(labels="Long sample name", axis=1)
    enzymeData = enzymeData[enzymeData["PlateCol"] <= 10]
    return enzymeData


def missingExtraMisnamed(enzymeData, dryWtSamples):
    """Checks the enzymeData dataframe to see if there are any plates that are
    missing, extra, or misnamed. Before running this function, the dataframe's
    Long sample names & Well columns must have been processed. After running
    this function and checking the output of this function, dry weight data
    can be added to the enzymeData dataframe.


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of enzyme activity that had undergone processing of its long
        sample name & well columns.
    dryWtSamples : List of strings
        List of strings where each string is a sample name taken from the dry
        weight spreadsheets. The names in the dry weight spreadsheets are
        correct; they are not misnamed, missing, or extra.

    Returns
    -------
    plateCounts : groupby() object
        This object holds the counts of plates to determine if there are any
        extra plates
    missingOrMisnamed : List of strings
        List of strings where each string is the name of a sample in
        dryWtSamples that is missing from the enzymeData dataframe. Can be used
        to check for any misnamed plates, too.
    """
    plateCounts = enzymeData.groupby("ID")["ID"].count()/96
    plateNames = plateCounts.index.tolist()
    samplesNotInPlateNames = []
    for sample in dryWtSamples:
        if sample not in plateNames:
            samplesNotInPlateNames.append(sample)
    return plateCounts, samplesNotInPlateNames


def addDryWtT035(enzymeData, processedDryWt, timepoint: str):
    """Adds dry weight data to the enzymeData dataframe. This function can only
    add the dry weight data to enzymeData dataframes from timepoints T0, T3,
    and T5. That is because the dry weight data for these 3 timepoints are from
    the same spreadsheet. The dry weight data for T6 had been processed in a
    different spreadsheet, so this function is not applicable to T6.

    This function can only be run after the enzymeData dataframe had been
    checked for any missing, misnamed, or extra plates. After this function had
    been ran, then treatment information can be added to the enzymeData
    dataframe.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of enzyme activity that had been checked for missing,
        misnamed, or extra plates
    processedDryWt : Pandas dataframe
        Dataframe of dry weight data for T0, T3, and T5 timepoints.
    timepoint : str
        String that describes the timepoint in enzymeData.

    Returns
    -------
    enzymeData : Pandas dataframe
        Dataframe of enzyme activity with dry weight data added.
    """
    timepointBool = processedDryWt["Time"].isin([timepoint])
    dryDF = processedDryWt[timepointBool]
    enzymeData = pd.merge(enzymeData, dryDF, how="left", on="ID")
    enzymeData = enzymeData.drop(labels="Time", axis=1)
    return enzymeData


def addTreatments(enzymeData):
    """Adds precipitation and vegetation treatments to the data frame of enzyme
    data based on sample names. Before this function can be called, the
    dataframe must have dry weight data.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of raw enzyme data where the "ID" column contains the sample
        ID. This column will be used to generate treatment information

    Returns
    -------
    The dataframe with 2 new columns (Precip & Vegetation) which represent
    treatments
    """
    for index, row in enzymeData.iterrows():
        rowID = row["ID"]
        if rowID != "B":
            precip = rowID[-2]
            plot = int(rowID[:-3])
            if plot <= 24:
                enzymeData.loc[index, "Vegetation"] = "Grassland"
            elif plot > 24:
                enzymeData.loc[index, "Vegetation"] = "CSS"
            if precip == "X":
                enzymeData.loc[index, "Precip"] = "Ambient"
            elif precip == "R":
                enzymeData.loc[index, "Precip"] = "Reduced"
            enzymeData.loc[index, "Plot"] = plot
    return enzymeData
# %%
# Following functions are reserved for black plates only. These functions are
# used to wrangle control data in black plates to make sure that control
# columns are side by side with the assay column


def subCtrlWrangling(enzymeData):
    """Manipulates buffer data so that they appear in a separate column, side
    by side with the Assay column in a black plate dataframe. Before this
    function can be called, the dataframe must have been processed so that its
    "Long sample name" and "Well" columns are each split into 2 new columns,
    and the "Long sample name" column is dropped. In addition, the dataframe
    must have dry weight data and treatment information.


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of hydrolytic enzyme activity that had undergone some
        processing

    Returns
    -------
    Dataframe with buffer readings now as a substrate control column.
    """
    bufferDF = enzymeData[enzymeData["ID"] == "B"]
    enzymeData = enzymeData[enzymeData["ID"] != "B"]
    colsToDrop = ["ID", "Dry assay (g)", "Vegetation", "Precip", "Plot"]
    bufferDF = bufferDF.drop(labels=colsToDrop, axis=1)
    bufferDF = bufferDF.rename(mapper={"Assay": "BufferReading"}, axis=1)
    subCtrlMerge = ["Well", "Assay date", "PlateRow", "PlateCol"]
    enzymeData = pd.merge(enzymeData, bufferDF, how="inner", on=subCtrlMerge)
    enzymeData = enzymeData.sort_values(by=["PlateCol", "Plot"])
    return enzymeData


def stanQuench(enzymeData):
    """Wrangles data in columns 8 and 9 so that they can serve as quench
    control and be used for their standard fluorescence. Columns 8 & 9 of each
    sample plate will be used as the sample plate's quench control, while the
    same columns of a buffer plate will be used as a sample plate's standard
    fluorescence.

    Before this function can be ran, the substrate control data in the
    enzymeData dataframe must have been manipulated to be in a separate column
    that is side by side with the "Assay" column.


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Hydrolytic enzyme activity dataframe with substrate control columns.

    Returns
    -------
    enzymeData : Pandas dataframe
        Hydrolytic enzyme activity dataframe with standard fluorescence &
        quench readings added to new columns.
    """
    AMC_DF = enzymeData[enzymeData["PlateCol"] == 8]
    AMC_DF = AMC_DF.drop(labels="Well", axis=1)
    MUB_DF = enzymeData[enzymeData["PlateCol"] == 9]
    MUB_DF = MUB_DF.drop(labels="Well", axis=1)
    homCtrlDF = enzymeData[enzymeData["PlateCol"] == 10]
    enzymeData = enzymeData[enzymeData["PlateCol"] <= 7]
    
    # Setting black plate columns that use a specific substrate for quench control
    MUB_DF1 = MUB_DF.copy()
    MUB_DF2 = MUB_DF.copy()
    MUB_DF3 = MUB_DF.copy()
    MUB_DF4 = MUB_DF.copy()
    MUB_DF5 = MUB_DF.copy()
    MUB_DF7 = MUB_DF.copy()
    
    standardFrames = [MUB_DF1, MUB_DF2, MUB_DF3, MUB_DF4, MUB_DF5, AMC_DF, MUB_DF7]
    for n in range(len(standardFrames)):
        currentDF = standardFrames[n]
        currentDF["PlateCol"] = n + 1
    
    # Merging AMC_DF & MUB_DF into a single dataframe of standard fluorescence
    # & quench control. This subsequent dataframe will be merged back into
    # enzymeData
    standardDF = pd.concat(objs=standardFrames, axis=0)
    standardDF = standardDF.sort_values(by=["PlateCol", "Plot"])
    # This is the right way to sort any dataframe that's derived from the enzyme
    # data to ensure that the sorted data frame looks like the original file
    
    # Renaming columns in the standard data frame and dropping redundant columns.
    # I deem certain columns as redundant because this dataframe will be merged
    # back into T0 Black, which already contains the information in the redundant
    # columns
    newColumnsDict = {"Assay": "QuenchCtrl",
                      "BufferReading": "StanFluo"}
    stanColsToDrop = ["Dry assay (g)", "Vegetation", "Precip"]
    standardDF = standardDF.rename(columns=newColumnsDict)
    standardDF = standardDF.drop(labels=stanColsToDrop, axis=1)
    
    # Merging standard dataframe back into enzymeData
    stanMergeLabels = ["Assay date", "ID", "PlateRow", "PlateCol", "Plot"]
    enzymeData = pd.merge(enzymeData, standardDF, how="inner", on=stanMergeLabels)
    enzymeData = enzymeData.sort_values(by=["Plot", "PlateCol"])  # sorting enzymeData
    enzymeData = enzymeData.rename(columns={"BufferReading": "SubCtrl"})
