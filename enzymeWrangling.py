# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 18:51:39 2020

@author: Brian Nhan Thien Chung
This is a module that is meant to facilitate the wrangling of enzyme activity
data and to help plot them.
"""
import pandas as pd
import matplotlib.pyplot as py
py.style.use("dark_background")
# %%
# Functions to wrangle both hydrolytic & oxidative enzyme data


def longNamesAndWells(enzymeData):
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


def dryWtT035(enzymeData, processedDryWt, timepoint: str):
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


def treatments(enzymeData):
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
    enzymeData = enzymeData.sort_values(by=["PlateCol", "Plot"])
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
    that is side by side with the "Assay" column. After this function had been
    ran, homogenate control data can be wrangled.


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Hydrolytic enzyme activity dataframe with substrate control columns.

    Returns
    -------
    enzymeData : Pandas dataframe
        Hydrolytic enzyme activity dataframe with standard fluorescence &
        quench readings added to new columns.
    homCtrlDF : Pandas dataframe
        Dataframe of homogenate control readings.
    """
    AMC_DF = enzymeData[enzymeData["PlateCol"] == 8]
    AMC_DF = AMC_DF.drop(labels="Well", axis=1)
    MUB_DF = enzymeData[enzymeData["PlateCol"] == 9]
    MUB_DF = MUB_DF.drop(labels="Well", axis=1)
    homCtrlDF = enzymeData[enzymeData["PlateCol"] == 10]
    enzymeData = enzymeData[enzymeData["PlateCol"] <= 7]

    # Setting black plate columns to their respective standards.
    MUB_DF1 = MUB_DF.copy()
    MUB_DF2 = MUB_DF.copy()
    MUB_DF3 = MUB_DF.copy()
    MUB_DF4 = MUB_DF.copy()
    MUB_DF5 = MUB_DF.copy()
    MUB_DF7 = MUB_DF.copy()

    standardFrames = [MUB_DF1, MUB_DF2, MUB_DF3,
                      MUB_DF4, MUB_DF5, AMC_DF, MUB_DF7]
    for n in range(len(standardFrames)):
        currentDF = standardFrames[n]
        currentDF["PlateCol"] = n + 1

    # Merging AMC_DF & MUB_DF into a single dataframe of standard fluorescence
    # & quench control. This subsequent dataframe will be merged back into
    # enzymeData
    standardDF = pd.concat(objs=standardFrames, axis=0)
    standardDF = standardDF.sort_values(by=["PlateCol", "Plot"])
    # This is the right way to sort any dataframe that's derived from the
    # enzyme data to ensure that the sorted data frame looks like the original
    # file

    # Renaming columns in the standard data frame and dropping redundant
    # columns. I deem certain columns as redundant because this dataframe will
    # be merged back into T0 Black, which already contains the information in
    # the redundant columns
    newColumnsDict = {"Assay": "QuenchCtrl",
                      "BufferReading": "StanFluo"}
    stanColsToDrop = ["Dry assay (g)", "Vegetation", "Precip"]
    standardDF = standardDF.rename(columns=newColumnsDict)
    standardDF = standardDF.drop(labels=stanColsToDrop, axis=1)

    # Merging standard dataframe back into enzymeData
    stanMerge = ["Assay date", "ID", "PlateRow", "PlateCol", "Plot"]
    enzymeData = pd.merge(enzymeData, standardDF, how="inner", on=stanMerge)
    enzymeData = enzymeData.sort_values(by=["Plot", "PlateCol"])
    enzymeData = enzymeData.rename(columns={"BufferReading": "SubCtrl"})
    return enzymeData, homCtrlDF


def homCtrlWrangling(enzymeData, homCtrlDF):
    """Wrangles homogenate control data. Column 10 of each sample plate will be
    used as that sample's homogenate control.

    Before this function can be run, standard fluorescence & quench control
    readings must have been wrangled. After this function is called, activity
    calculations can take place.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of enzyme activity data after standard fluorescence & quench
        control readings have been wrangled.
    homCtrlDF : Pandas dataframe
        Homogenate control readings.

    Returns
    -------
    enzymeData : Pandas dataframe
        Dataframe of enzyme activity data with homogenate control data added
    """
    homCtrlDF = homCtrlDF.rename(mapper={"Assay": "HomCtrl"}, axis=1)
    homCtrlColsToDrop = ["Well", "Assay date", "PlateCol", "Dry assay (g)",
                         "Vegetation", "Precip", "BufferReading"]
    homCtrlDF = homCtrlDF.drop(labels=homCtrlColsToDrop, axis=1)
    homCtrlMerge = ["ID", "Plot", "PlateRow"]
    enzymeData = pd.merge(enzymeData, homCtrlDF, how="inner", on=homCtrlMerge)
    enzymeData = enzymeData.sort_values(by=["Plot", "PlateCol"])
    return enzymeData


def hydrolyticEnzymeActivity(enzymeData, plateInfo):
    """Calculates hydrolytic enzyme activity using formulas from German et al
    2011.

    This function can only be run after all of the control readings have been
    properly wrangled so that they appear in separate columns side by side with
    the Assay column. Once this is ran, then the enzyme activities in the
    subsequent dataframe can be plotted.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Pandas dataframe with all the control readings wrangled.
    plateInfo : Pandas dataframe
        Pandas dataframe that contains miscellaneous plate information: amount
        of standards in each well in columns 8 & 9, substrate concentrations,
        and enzyme names.

    Returns
    -------
    enzymeData : Pandas dataframe
        Pandas dataframe with calculated enzyme activities as well as added
        information from plateInfo.
    """
    # Calculating Quench Coefficients. Each quench coefficient is specific
    # to a particular row of a sample plate. Quench Coefficients are unitless
    enzymeData["QuenchCoef"] = ((enzymeData["QuenchCtrl"]
                                 - enzymeData["HomCtrl"])
                                / enzymeData["StanFluo"])

    # Calculating Emission Coefficients. Each coefficient is specific to a row
    # of a plate. Units: nmol^-1
    assayVol = 0.250  # mL. Volume of assay well, which consists of
    # substrate (0.125 mL) + homogenate (0.125 mL)
    stanVol = 0.250  # mL. Volume of quench and standard wells, which
    # consists of standard (0.125 mL) +
    # either homogenate(sample plate, 0.125 mL) or
    # buffer(buffer plate, 0.125 mL)
    enzymeData["EmisCoef"] = ((enzymeData["StanFluo"]*stanVol)
                              / (enzymeData["StanAmt"]*assayVol))

    # Calculating net fluorescence. Units are fluorescence units, which is
    # essentially meaningless and so can be considered as unitless
    enzymeData["NetFluo"] = (((enzymeData["Assay"] - enzymeData["HomCtrl"])
                              / enzymeData["QuenchCoef"])
                             - enzymeData["SubCtrl"])

    # Calculating enzyme activity. Initial units are nmol L^-1 g^-1 h^-1. Final
    # units will be micromole L^-1 g^-1 h^-1. This enzyme activity had been
    # normalized for the mass of litter used in the assay
    incubaTime = 4  # Hours. Each plate was left to sit for 4 hours
    bufferVol = 150  # mL. Volume of buffer used in preparing a homogenate
    homVol = 0.125  # mL. Volume of homogenate that is pipetted into each well
    # in columns 1-7 in each sample plate.
    enzymeData["Activity"] = ((enzymeData["NetFluo"]*bufferVol)
                              / (enzymeData["EmisCoef"]*homVol*incubaTime
                              * enzymeData["Dry assay (g)"]))
    enzymeData["Activity"] = enzymeData["Activity"]/1000
    enzymeData = enzymeData.sort_values(by=["Plot", "PlateCol"])
    return enzymeData


def plotHydrolyticEnzymeActivity(enzymeData, dryWtSamples, plotPath, tpoint):
    """Plots the hydrolytic enzyme activity of a dataframe that represents
    a timepoint and contains calculated hydrolytic enzyme activity in that
    timepoint

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of calculated hydrolytic enzyme activity from a single
        timepoint
    dryWtSamples : List of strings
        List of strings, with each string being the name of a sample from a
        timepoint. I use the sample names generated from the dry weight
        spreadsheet for T0, T3, and T5 because I know that these names are not
        missing, misnamed, or extras
    plotPath : Path object
        The path in which the plots of hydrolytic enzyme activity will be saved
    tpoint : str
        String representing the timepoint of the enzymeData dataframe.

    Returns
    -------
    None.
    """
    enzymes = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG"]
    for sample in dryWtSamples:
        sampleDF = enzymeData[enzymeData["ID"] == sample]
        vegetation = sampleDF.groupby("Vegetation")["Vegetation"].count().index
        vegetation = vegetation.tolist()
        vegetation = vegetation[0]
        precip = sampleDF.groupby("Precip")["Precip"].count().index.tolist()
        precip = precip[0]
        figTitle = "{0:}, {1:}, {2:}, {3:}".format(sample, vegetation, precip,
                                                   tpoint)
        py.figure(num=figTitle, figsize=(25, 10))
        for enzyme in enzymes:
            if enzyme == "NAG":
                plotIndex = 8
            else:
                plotIndex = enzymes.index(enzyme) + 1
            py.subplot(2, 4, plotIndex)
            substrateDF = sampleDF[sampleDF["Enzyme"] == enzyme]
            py.plot(substrateDF["SubConcen"], substrateDF["Activity"])
            py.scatter(substrateDF["SubConcen"], substrateDF["Activity"])
            py.title("{0:}, {1:}, {2:}, {3:}".format(sample, vegetation,
                                                     precip, enzyme))
            py.xlabel("Substrate concentration (micromolar)")
            py.ylabel("Activity (micromole L^-1 g^-1 h^-1)")
        # Saving figures for data quality control purposes
        figName = figTitle + ".png"
        figPath = plotPath/figName
        py.savefig(figPath)
    return
