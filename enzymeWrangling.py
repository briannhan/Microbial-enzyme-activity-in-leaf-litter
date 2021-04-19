# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 18:51:39 2020

@author: Brian Nhan Thien Chung
This is a module that is meant to facilitate the wrangling of enzyme activity
data and to help plot them.
"""
import pandas as pd
import matplotlib.pyplot as py
import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
py.style.use("dark_background")

# %%
# Making dataframes of substrate concentrations, enzyme names, amounts of
# standards, and other plate information


# Plate information for black plates
AMCamount = 62.5*125/1000
"""Amount of MUB standard in standard and quench control wells.
62.5 micromolar (concentration) x 125 microliter pipetted/ 1000 (conversion)
units: nanomoles"""
MUBamount = 25*125/1000
"""Amount of AMC standard in standard and quench control wells.
25 micromolar (concentration) x 125 microliter pipetted/ 1000 (conversion)
units: nanomoles"""
plateCols = np.linspace(start=1, stop=7, num=7).tolist()
enzymeName = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG"]
standardAmount = [MUBamount, MUBamount, MUBamount, MUBamount, MUBamount,
                  AMCamount, MUBamount]
hydroInfoDic = {"PlateCol": plateCols, "Enzyme": enzymeName,
                "StanAmt": standardAmount}
hydroInfo = pd.DataFrame(hydroInfoDic)
# Making dataframe of hydrolytic enzyme concentrations to merge back into
# hydroInfo.
AGnames = 8*["AG"]
APnames = 8*["AP"]
BGnames = 8*["BG"]
BXnames = 8*["BX"]
CBHnames = 8*["CBH"]
LAPnames = 8*["LAP"]
NAGnames = 8*["NAG"]
enzymesLists = [AGnames, APnames, BGnames, BXnames,
                CBHnames, LAPnames, NAGnames]
longEnzymesNames = [name for outerList in enzymesLists for name in outerList]
subProps = np.geomspace(1, 1/128, num=8)
# Units of substrate concentrations are in micromolar
AGconcen = (500*subProps).tolist()
APconcen = (2000*subProps).tolist()
BGconcen = (500*subProps).tolist()
BXconcen = (500*subProps).tolist()
CBHconcen = (250*subProps).tolist()
LAPconcen = (500*subProps).tolist()
NAGconcen = (1000*subProps).tolist()
subConcenLists = [AGconcen, APconcen, BGconcen, BXconcen, CBHconcen,
                  LAPconcen, NAGconcen]
subConcen = [concen for outer in subConcenLists for concen in outer]
plateRow = 7*list("ABCDEFGH")
blackSubConcenDict = {"PlateRow": plateRow, "Enzyme": longEnzymesNames,
                      "SubConcen": subConcen}
blackSubConcenDF = pd.DataFrame(blackSubConcenDict)
hydroInfo = pd.merge(hydroInfo, blackSubConcenDF, how="inner", on="Enzyme")


# Plate information for clear plates
pyroHighest = (1e6)/(7.9*126.11*2)
'''Highest concentration of pyrogallol is 1 mg pyrogallol/7.9 mL water. Molar
mass of pyrogallol is 126.11 g/mol. I multiplied by 1,000,000 to give the final
value units of micromole L^-1. I divided by 2 to take into account the
fact that half the volume of each assay well consists of the pipetted
substrate (pyrogallol) and the other half consists of the filtered
homogenate.'''
pyroConcen = (pyroHighest*subProps).tolist()
pyroConcenDF = pd.DataFrame({"PlateRow": list("ABCDEFGH"),
                             "SubConcen": pyroConcen})
# %%
# Functions to wrangle both hydrolytic & oxidative enzyme data


def longNamesAndWells(enzData):
    """Processes the long sample names & wells in an enzyme dataframe that
    occurred from reading in a dataframe of raw enzyme data. The long sample
    names are in a column, which are then split up into 2 columns. The Well
    column will also be used to create 2 new columns


    Parameters
    ----------
    enzData : Pandas dataframe
        Dataframe of raw enzyme data from reading in a raw enzyme data text
        file. Has 3 columns: Well, Long sample name, Assay. Well will be used
        to create 2 new columns, while Long sample name will also be used
        to create 2 new columns

    Returns
    -------
    The dataframe with the names & wells processed & with new columns. Now,
    the dataframe can be checked for missing, extra, misnamed plates

    """
    enzData["Long sample name"] = enzData["Long sample name"].str.split("_")
    for index, row in enzData.iterrows():
        longSampleName = row["Long sample name"]  # splitted list
        enzData.loc[index, "Assay date"] = longSampleName[0]
        for part in longSampleName:
            if "X" in part:
                sampleID = part
        if "B" in longSampleName:
            enzData.loc[index, "ID"] = "B"
        elif "X" in sampleID:
            enzData.loc[index, "ID"] = sampleID
        well = row["Well"]
        wellRow = well[0]
        wellColumn = int(well[1:])
        enzData.loc[index, "PlateRow"] = wellRow
        enzData.loc[index, "PlateCol"] = wellColumn
    enzData = enzData.drop(labels="Long sample name", axis=1)
    return enzData


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


def dryWt(enzymeData, processedDryWt, timepoint: str):
    """Adds dry weight data to the enzymeData dataframe. While the dry weight
    for T6 is originally recorded in a different spreadsheet than the
    spreadsheet that contains dry weight for T0, T3, & T5, dry weight from all
    timepoints are combined into a single dataframe, so this function is
    applicable to all timepoints.

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
    enzymeData = enzymeData.dropna(axis=0, how="all")
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
        if rowID != "B" and type(rowID) == str:
            # Added type of rowID conditions because of error in processing
            # T5 black plates. This error is due to a difference in naming
            # schemes between some samples in T5 and previous timepoints T0 &
            # T3
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
        elif type(rowID) != str:
            print(type(rowID))
            print(row)
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


def hydrolaseActivity(enzymeData):
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

    Returns
    -------
    enzymeData : Pandas dataframe
        Pandas dataframe with calculated enzyme activities as well as added
        information from hydroInfo.
    """
    # Calculating Quench Coefficients. Each quench coefficient is specific
    # to a particular row of a sample plate. Quench Coefficients are unitless
    enzymeData = pd.merge(enzymeData, hydroInfo, on=["PlateCol", "PlateRow"])
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


def plotHydrolaseActivity(enzymeData, plotPath, tpoint):
    """Plots the hydrolytic enzyme activity of a dataframe that represents
    a timepoint and contains calculated hydrolytic enzyme activity in that
    timepoint

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of calculated hydrolytic enzyme activity from a single
        timepoint
    plotPath : Path object
        The path in which the plots of hydrolytic enzyme activity will be saved
    tpoint : str
        String representing the timepoint of the enzymeData dataframe.

    Returns
    -------
    None.
    """
    enzymes = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG"]
    samples = enzymeData.groupby("ID")["ID"].count().index.tolist()
    for sample in samples:
        sampleDF = enzymeData[enzymeData["ID"] == sample]
        vegetation = sampleDF.groupby("Vegetation")["Vegetation"].count().index
        vegetation = vegetation.tolist()
        vegetation = vegetation[0]
        precip = sampleDF.groupby("Precip")["Precip"].count().index.tolist()
        precip = precip[0]
        date = sampleDF.groupby("Assay date")["Assay date"].count().index
        date = date.tolist()
        date = date[0]
        figTitle = "{0}, {1}, {2}, {3}, {4}".format(sample, vegetation, precip,
                                                    tpoint, date)
        py.figure(num=figTitle, figsize=(25, 10))
        for enzyme in enzymes:
            if enzyme == "NAG":
                plotIndex = 8
            else:
                plotIndex = enzymes.index(enzyme) + 1
            py.subplot(2, 4, plotIndex)
            substrateDF = sampleDF[sampleDF["Enzyme"] == enzyme]
            py.plot("SubConcen", "Activity", data=substrateDF)
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
# %%
# This section wrangles and plots enzyme activity of clear plate data


def enzymesNreplicates(enzymeData):
    """Assigns enzyme names and replicates to the enzymeData dataframe. This
    function must be ran first before the dataframe can be split into
    dataframes of buffer readings, substrate and homogenate controls, and
    assay readings.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        A dataframe of oxidative enzyme data.

    Returns
    -------
    enzymeData : Pandas dataframe
        The dataframe with the enzyme names and replicate numbers added.
    """
    PPOcols = [1, 4, 5, 8, 9, 12]
    PERPPOcols = [2, 3, 6, 7, 10, 11]
    replicate1 = [5, 6, 7, 8]
    replicate2 = [9, 10, 11, 12]
    for index, row in enzymeData.iterrows():
        if row["PlateCol"] in PPOcols:
            enzymeData.loc[index, "Enzyme"] = "PPO"
        elif row["PlateCol"] in PERPPOcols:
            enzymeData.loc[index, "Enzyme"] = "Both"
        if row["PlateCol"] in replicate1:
            enzymeData.loc[index, "Replicate"] = 1
        elif row["PlateCol"] in replicate2:
            enzymeData.loc[index, "Replicate"] = 2
    return enzymeData


def separateControls(enzymeData):
    """Separates the enzymeData dataframe into separate dataframes comprising
    of absorbance readings of the buffer, the substrate and homogenate
    controls, and the assay readings.

    Before this function can be run, the enzymeData dataframe must have 2
    columns that dictate the enzyme and replicate. After this function is run,
    the control readings and assay readings can be blanked.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        A dataframe of oxidative enzyme data.

    Returns
    -------
    bufferAbs : Pandas dataframe
        A dataframe of the absorbance of the buffer wells (which contains
        buffer and water)
    subConAbs : Pandas dataframe
        A dataframe of the absorbance of the substrate control wells (which
        contains buffer and substrate)
    homCtrlAbs : Pandas dataframe
        A dataframe of the absorbance of the homogenate control wells (which
        contains the filtered homogenate and water)
    enzymeData : Pandas dataframe
        The original input dataframe with the buffer and substrate and
        homogenate control readings removed
    """
    bufferBool = enzymeData["PlateCol"].isin([3, 4])
    bufferAbs = enzymeData[bufferBool]
    subConAbsBool = enzymeData["PlateCol"].isin([1, 2])
    subConAbs = enzymeData[subConAbsBool]
    homCtrlAbsBool = enzymeData["PlateCol"].isin([7, 8, 11, 12])
    homCtrlAbs = enzymeData[homCtrlAbsBool]
    assayAbsBool = enzymeData["PlateCol"].isin([5, 6, 9, 10])
    enzymeData = enzymeData[assayAbsBool]

    ctrlColsToDrop = ["PlateCol", "Dry assay (g)", "Vegetation", "Precip"]
    bufferAbs = bufferAbs.drop(labels=ctrlColsToDrop, axis=1)
    subConAbs = subConAbs.drop(labels=ctrlColsToDrop, axis=1)
    homCtrlAbs = homCtrlAbs.drop(labels=ctrlColsToDrop, axis=1)
    return bufferAbs, subConAbs, homCtrlAbs, enzymeData


def blank(enzymeData):
    """Removes the absorbance of the buffer from the substrate and homogenate
    controls and the assay wells. I'm calling this process 'blanking'.

    This function calls the separateControls() function that separate control,
    buffer, and assay readings readings into separate dataframes. After this
    function is run, the control readings can be merged back into enzymeData.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        A dataframe of oxidative enzyme activity data

    Returns
    -------
    subConAbs : Pandas dataframe
        A dataframe of the absorbance of the substrate control wells (which
        contains buffer and substrate) which the absorbance of the buffer
        removed
    homCtrlAbs : Pandas dataframe
        A dataframe of the absorbance of the homogenate control wells (which
        contains the filtered homogenate and water) with the absorbance of the
        buffer removed
    enzymeData : Pandas dataframe
        A dataframe with only assay readings with the absorbance of the buffer
        removed
    """
    bufferAbs, subConAbs, homCtrlAbs, enzymeData = separateControls(enzymeData)
    bufferAbs = bufferAbs.rename(mapper={"Assay": "Buffer"}, axis=1)
    bufferAbsColsToDrop = ["Well", "Assay date", "Replicate"]
    bufferAbs = bufferAbs.drop(labels=bufferAbsColsToDrop, axis=1)
    bufferAbsMerge = ["ID", "Plot", "PlateRow", "Enzyme"]
    subConAbs = pd.merge(left=subConAbs, right=bufferAbs,
                         how="inner", on=bufferAbsMerge)
    homCtrlAbs = pd.merge(left=homCtrlAbs, right=bufferAbs,
                          how="inner", on=bufferAbsMerge)
    enzymeData = pd.merge(left=enzymeData, right=bufferAbs, how="inner",
                          on=bufferAbsMerge)
    subConAbs["Assay"] = subConAbs["Assay"] - subConAbs["Buffer"]
    homCtrlAbs["Assay"] = homCtrlAbs["Assay"] - homCtrlAbs["Buffer"]
    enzymeData["Assay"] = enzymeData["Assay"] - enzymeData["Buffer"]
    return subConAbs, homCtrlAbs, enzymeData


def processCtrls(enzymeData):
    """Merge control dataframes into the enzymeData dataframe to produce a
    single dataframe with assay readings and control readings side by side to
    facilitate calculations.

    This function calls the blank() function that removes the absorbance of
    buffer from the control readings and the assay readings. After this
    function is run, oxidase activities can be calculated.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        A dataframe of oxidative enzyme absorbance readings.

    Returns
    -------
    enzymeData : Pandas dataframe
        The original enzymeData dataframe with substrate and homogenate control
        readings added, so that this dataframe has assay readings and readings
        from substrate and homogenate controls and all readings are blanked.
    """
    subConAbs, homCtrlAbs, enzymeData = blank(enzymeData)
    subConAbsColsToDrop = ["Well", "Assay date", "Replicate", "Buffer"]
    subConAbs = subConAbs.rename(columns={"Assay": "SubCtrl"})
    subConAbs = subConAbs.drop(labels=subConAbsColsToDrop, axis=1)
    homCtrlAbsColsToDrop = ["Well", "Assay date", "Buffer"]
    homCtrlAbs = homCtrlAbs.rename(columns={"Assay": "HomCtrl"})
    homCtrlAbs = homCtrlAbs.drop(labels=homCtrlAbsColsToDrop, axis=1)
    subConAbsMerge = ["ID", "Plot", "PlateRow", "Enzyme"]
    enzymeData = pd.merge(left=enzymeData, right=subConAbs, how="inner",
                          on=subConAbsMerge)
    homCtrlAbsMerge = ["ID", "PlateRow", "Plot", "Enzyme", "Replicate"]
    enzymeData = pd.merge(left=enzymeData, right=homCtrlAbs, how="inner",
                          on=homCtrlAbsMerge)
    enzymeData = enzymeData.drop(labels="Buffer", axis=1)
    enzymeData = enzymeData.sort_values(by=["PlateCol", "Plot"])
    return enzymeData


def oxidaseActivity(enzymeData):
    """Calculates the oxidative enzyme activity of polyphenol oxidase (PPO)
    and peroxidase (PER)

    Before this function can be run, the input dataframe must contain blanked
    out substrate and homogenate control readings. After this function is run,
    then the oxidative enzyme activities can be plotted.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of oxidative enzyme data, containing blanked out assay
        readings and blanked out control readings.

    Returns
    -------
    enzymeData : Pandas dataframe
        Dataframe of oxidative enzyme data that now contain PPO and PER
        activities
    """
    extincCoef = 4.2  # micromole^-1
    bufferVol = 150  # mL. Volume of buffer to prepare a single homogenate
    homVol = 0.125  # mL. Volume of homogenate that is pipetted into each well
    enzymeData["NetAbs"] = (enzymeData["Assay"] - enzymeData["HomCtrl"]
                            - enzymeData["SubCtrl"])
    incubaTimeClear = 24  # hours
    enzymeData["Activity"] = ((enzymeData["NetAbs"]*bufferVol) /
                              (extincCoef*homVol*incubaTimeClear
                               * enzymeData["Dry assay (g)"]))
    # Units of enzyme activity: micromole g^-1 hr^-1

    # Wrangle calculated enzyme activities to extract PER from combined PER &
    # PPO activities
    allOxi = enzymeData[enzymeData["Enzyme"] == "Both"]
    enzymeData = enzymeData[enzymeData["Enzyme"] == "PPO"]
    allOxi = allOxi.rename(columns={"Activity": "All oxidases"})
    enzymeData = enzymeData.rename(columns={"Activity": "PPO activity"})
    allOxiLabelsToDrop = ["Well", "Assay", "Assay date", "PlateCol",
                          "Dry assay (g)", "Vegetation", "Precip", "Enzyme",
                          "SubCtrl", "HomCtrl", "NetAbs"]
    allOxi = allOxi.drop(columns=allOxiLabelsToDrop)
    enzymeData = enzymeData.drop(columns="Enzyme")
    allOxiMerge = ["ID", "PlateRow", "Plot", "Plot", "Replicate"]
    enzymeData = pd.merge(left=enzymeData, right=allOxi, how="inner",
                          on=allOxiMerge)

    # Extract PER from combined PER & PPO activities
    enzymeData = enzymeData.rename(columns={"All oxidases": "PER activity"})
    enzymeData["PER activity"] = (enzymeData["PER activity"]
                                  - enzymeData["PPO activity"])
    enzymeData = enzymeData.sort_values(by=["PlateCol", "Plot"])
    return enzymeData


def plotOxidaseActivity(enzymeData, plotPath):
    """Plots oxidase activities.

    Before this function can be run, the enzymeData dataframe must contain
    calculated enzyme activities.

    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of calculated enzyme activity
    plotPath : Path object
        The path in which the plots of oxidase activity will be saved

    Returns
    -------
    None.
    """
    enzymeData = pd.merge(left=enzymeData, right=pyroConcenDF, how="inner",
                          on="PlateRow")
    samples = enzymeData.groupby("ID")["ID"].count().index.tolist()
    # Graph oxidative enzyme activity
    for sample in samples:
        sampleDF = enzymeData[enzymeData["ID"] == sample]
        vegetation = sampleDF.groupby("Vegetation")["Vegetation"].count().index
        vegetation = vegetation.tolist()
        vegetation = vegetation[0]
        precip = sampleDF.groupby("Precip")["Precip"].count().index.tolist()
        precip = precip[0]
        date = sampleDF.groupby("Assay date")["Assay date"].count().index
        date = date.tolist()
        date = date[0]
        figureTitle = "{0}, {1}, {2}, {3}".format(sample, vegetation,
                                                  precip, date)
        py.figure(num=figureTitle, figsize=(17, 13))
        for i in range(2):
            # Making subplot of PPO for 1 replicate
            replicate = i + 1
            py.subplot(2, 2, (i*2) + 1)
            replicateDF = sampleDF[sampleDF["Replicate"] == replicate]
            py.title("{0}, {1}, {2}, PPO, replicate {3}".format(sample,
                                                                vegetation,
                                                                precip,
                                                                replicate))
            py.xlabel("Substrate concentration (micromole L^-1)")
            py.ylabel("Normalized enzyme activity (micromole g^-1 hr^-1)")
            py.plot("SubConcen", "PPO activity", data=replicateDF)
            py.scatter("SubConcen", "PPO activity", data=replicateDF)

            # Making subplot of PER for 1 replicate
            py.subplot(2, 2, (i*2) + 2)
            py.title("{0}, {1}, {2}, PER, replicate {3}".format(sample,
                                                                vegetation,
                                                                precip,
                                                                replicate))
            py.xlabel("Substrate concentration (micromole L^-1)")
            py.ylabel("Normalized enzyme activity (micromole g^-1 hr^-1)")
            py.plot("SubConcen", "PER activity", data=replicateDF)
            py.scatter("SubConcen", "PER activity", data=replicateDF)
        figureName = figureTitle + ".png"
        figPath = plotPath/figureName
        py.savefig(figPath)
    return


# %%
# Cleaning calculated enzyme activity by (1) setting negative activity values
# to 0 and (2) removing substrate inhibition


def cleanKeysDFs(processDF):
    """
    Produces 2 dataframes -- 1 for hydrolase activity, 1 for oxidase
    activity -- of cleaning keys that will be used to clean dataframes of
    calculated hydrolase or oxidase activity

    Parameters
    ----------
    processDF : Pandas dataframe
        Original dataframe of processing keys as read in. This dataframe will
        be separated into 2 dataframes.

    Returns
    -------
    2 dataframes of keys, 1 for hydrolase activity, 1 for oxidase activity

    """
    columns = processDF.columns.tolist()
    # The last column is called "Unnamed" as in the original Excel file, it
    # does not have a name. This column will be removed
    columns = columns[:-1]

    hydroEnzymes = columns[0:7]
    oxidaseEnzymesNreplicates = columns[7:]
    hydroCleanDF = processDF.copy()
    oxiCleanDF = processDF.copy()  # haha
    hydroCleanDF = hydroCleanDF.drop(columns=oxidaseEnzymesNreplicates)
    oxiCleanDF = oxiCleanDF.drop(columns=hydroEnzymes)
    return hydroCleanDF, oxiCleanDF


def cleanHydro(processDF, data):
    """
    Processes a hydrolytic enzyme activity dataframe by setting negative
    activity values to 0 and removing substrate inhibition. These are the 2
    main problems that plague the data quality of the calculated hydrolase
    activity, although relatively few samples have the problem of negative
    hydrolase activity.

    Parameters
    ----------
    processDF : Pandas dataframe
        Dataframe containing the keys that will be used to denote the types
        of processing to be applied to a sample's enzyme. This dataframe will
        be the original dataframe as read in from the "Samples with errors"
        Excel file, and will be processed using the cleanKeysDFs() function.
    data : Pandas dataframe
        Dataframe of hydrolytic enzyme activity data. This dataframe had
        already been sorted so that for each enzyme of a sample, the rows are
        ordered in descending order by substrate concentration so that as you
        go down the rows of a particular sample's enzyme, the substrate
        concentration decreases.

    Returns
    -------
    The cleaned dataframe of hydrolytic enzyme activity data.

    """
    hydroProcessDF, oxiProcessDF = cleanKeysDFs(processDF)

    # Setting negative activity values to 0
    hydroProcessInd = data.index[data["Activity"] < 0].tolist()
    data.loc[hydroProcessInd, "Activity"] = 0

    # Removing substrate inhibition
    hydroEnzymes = hydroProcessDF.columns.tolist()[1:]
    initialShape = data.shape
    print("Initial shape is", initialShape)
    for index, row in hydroProcessDF.iterrows():
        for enzyme in hydroEnzymes:
            if row[enzyme] == "o":  # "o" stands for substrate inhibition
                sampleDF = data[data["ID"] == row["ID"]]
                dfToProcess = sampleDF[sampleDF["Enzyme"] == enzyme]
                indexMaxActivity = dfToProcess["Activity"].idxmax()
                dfIndices = dfToProcess.index.tolist()
                indicesToDrop = np.arange(dfIndices[0], indexMaxActivity)
                data = data.drop(index=indicesToDrop)
    finalShape = data.shape
    print("Final shape is", finalShape)
    return data
# %%
# Nonlinear regressions


def MM(S, Vmax, Km):
    """
    A Michaelis-Menten function. The function outputs the activity of a
    particular enzyme given a set of non-normalized substrate concentrations
    and the Michaelis-Menten parameters for that enzyme. The purpose of this
    function is to be fitted to enzyme activity data using nonlinear regression
    to produce parameter values that will then be statistically analyzed.

    Parameters
    ----------
    S : float/integers or Pandas series of floats/integers
        The non-normalized substrate concentration in units of micromolar. The
        normalized concentration would be if this is divided by dry litter
        mass.
    Vmax : float
        Normalized maximum reaction velocity, normalized because if nonlinear
        regression is used to fit this function, it will produce a Vmax that is
        divided by the litter mass. Units are micromole g^-1 hr^-1.
    Km : float
        Michaelis-Menten constant. While normally referred to as the
        half-saturation constant, the concentration of substrates at which
        the reaction velocity is half that of the maximum reaction velocity,
        its more precise definition is the equilibrium between the breakdown
        of the enzyme-substrate complex into products and substrates and the
        formation of the enzyme-substrate complex from substrates. It is the
        ratio between these rates of breakdown and formation. Units are
        micromolar, like substrate concentrations. The seemingly weird units
        of Vmax does not actually require any unit conversions for Km, so we
        can keep the same units on Km.

    Returns
    -------
    Reaction velocity, enzyme activity, units of micromole g^-1 hr^-1, similar
    to Vmax.

    """
    return (Vmax*S)/(Km + S)


def nonlinRegress(data, enzymeType, timepoint=None):
    """
    Performs nonlinear regression by fitting the MM function to a dataframe
    that contains the input for MM (substrate concentration) and the output
    values for MM (reaction velocity). As a dataframe of hydrolytic enzymes
    vs oxidative enzymes are formatted differently from each other, there will
    be some slight differences in how they are handled. Ultimately, this
    function produces a dataframe containing the enzyme parameters as output.

    Parameters
    ----------
    data : Pandas dataframe
        Dataframe of enzyme activity and substrate concentration.
    enzymeType : String
        Specifies the enzyme type (hydrolotic or oxidative) as there will be
        differences in how each are handled. Use 'H' for hydrolytic enzymes,
        aka black plates, and 'O' for oxidative enzymes, aka clear plates.
    timepoint : int
        The timepoint of when the data was collected (either 0, 3, 5, or 6).
        The purpose of this parameter is only to allow for the processing of
        sample 47RRX in timepoint 5, which was assayed twice for oxidative
        enzymes.

    Returns
    -------
    Pandas dataframe of enzyme parameters
    """
    # Obtaining hydrolytic enzyme names if data is hydrolase activity
    if enzymeType == "H":
        hydroEnzymes = data.groupby("Enzyme")["Enzyme"].count().index.tolist()
    elif enzymeType == "O" and timepoint == 5:
        # removing sample 47RRX because in T5, its oxidative enzyme activity
        # was assayed twice.
        samp47 = data[data["ID"] == "47RRX"]
        samp47_190125 = samp47[samp47["Assay date"] == "190125"]
        samp47_190222 = samp47[samp47["Assay date"] == "190222"]
        data = data[data["ID"] != "47RRX"]
    samples = data.groupby("ID")["ID"].count().index.tolist()

    VmaxDF = pd.DataFrame({"ID": samples})
    KmDF = pd.DataFrame({"ID": samples})

    # Performing regression for oxidase activity (if it's T5, then it doesn't
    # include 47RRX) as well as hydrolase activity
    for sample in samples:
        sampleDF = data[data["ID"] == sample]
        sampleIndex = VmaxDF[VmaxDF["ID"] == sample].index

        # Performing nonlinear regression for hydrolase activity
        if enzymeType == "H":
            for enzyme in hydroEnzymes:
                enzymeDF = sampleDF[sampleDF["Enzyme"] == enzyme]
                enzymeDF = enzymeDF.sort_values(by="SubConcen")
                try:
                    params, paramVar = curve_fit(MM, enzymeDF["SubConcen"],
                                                 enzymeDF["Activity"],
                                                 bounds=(0, np.inf))
                    VmaxDF.loc[sampleIndex, enzyme] = params[0]
                    KmDF.loc[sampleIndex, enzyme] = params[1]
                except RuntimeError:
                    print("Can't fit {0:}, enzyme {1:}".format(sample, enzyme))
                    VmaxDF.loc[sampleIndex, enzyme] = "can't fit"
                    KmDF.loc[sampleIndex, enzyme] = "can't fit"

        # Performing nonlinear regression for oxidase activity (not including
        # T5 sample 47RRX)
        elif enzymeType == "O":
            for i in range(2):
                replicate = i + 1
                repDF = sampleDF[sampleDF["Replicate"] == replicate]
                repDF = repDF.sort_values(by="SubConcen")
                try:  # fitting MM to PPO
                    paramsPPO, paramVarPPO = curve_fit(MM, repDF["SubConcen"],
                                                       repDF["PPO activity"],
                                                       bounds=(0, np.inf))
                    PPOcol = "PPO {0:d}".format(replicate)
                    VmaxDF.loc[sampleIndex, PPOcol] = paramsPPO[0]
                    KmDF.loc[sampleIndex, PPOcol] = paramsPPO[1]
                except RuntimeError:
                    print("Can't fit {1}, PPO replicate {0}".format(replicate,
                                                                    sample))
                    VmaxDF.loc[sampleIndex, PPOcol] = "can't fit"
                    KmDF.loc[sampleIndex, PPOcol] = "can't fit"

                try:  # fitting MM to PER
                    paramsPER, paramVarPER = curve_fit(MM, repDF["SubConcen"],
                                                       repDF["PER activity"],
                                                       bounds=(0, np.inf))
                    PERcol = "PER {0:d}".format(replicate)
                    VmaxDF.loc[sampleIndex, PERcol] = paramsPER[0]
                    KmDF.loc[sampleIndex, PERcol] = paramVarPER[1]
                except RuntimeError:
                    print("Can't fit {1}, PER replicate {0}".format(replicate,
                                                                    sample))
                    VmaxDF.loc[sampleIndex, PERcol] = "can't fit"
                    KmDF.loc[sampleIndex, PERcol] = "can't fit"

    # Performing nonlinear regression for T5 oxidase sample 47RRX
    if enzymeType == "O" and timepoint == 5:
        samp47frames = [samp47_190125, samp47_190222]
        sampleIndex = VmaxDF[VmaxDF["ID"] == "47RRX"].index
        for sampleDF in samp47frames:
            # Following for loop loops through the replicates for sample 47RRX.
            # The way that the parameters dataframes and the calculated enzyme
            # activities dataframe record replicates is different. In the
            # calculated enzyme activities dataframe, replicates are
            # either 1 or 2 even if this sample was assayed twice. In the
            # parameters dataframes, I intend the replicates to be a range of
            # integers with the following mathematical notation [1, 4].
            for i in range(2):
                replicateData = i + 1  # Replicate recorded on the dataframe
                # of calculated enzyme activities

                # Following if statements calculates replicate for the
                # parameter dataframes
                if sampleDF == samp47_190125:
                    repPar = i + 1
                elif sampleDF == samp47_190222:
                    repPar = i + 3

                repDF = sampleDF[sampleDF["Replicate"] == replicateData]
                repDF = repDF.sort_values(by="SubConcen")
                try:  # Fitting MM to PPO activity
                    paramsPPO, paramVarPPO = curve_fit(MM, repDF["SubConcen"],
                                                       repDF["PPO activity"],
                                                       bounds=(0, np.inf))
                    PPOcol = "PPO {0:d}".format(repPar)
                    VmaxDF.loc[sampleIndex, PPOcol] = paramsPPO[0]
                    KmDF.loc[sampleIndex, PPOcol] = paramsPPO[1]
                except RuntimeError:
                    print("Can't fit {0:}, PPO replicate {1:}".format(sample,
                                                                      repPar))
                    VmaxDF.loc[sampleIndex, PPOcol] = "can't fit"
                    KmDF.loc[sampleIndex, PPOcol] = "can't fit"

                try:  # fitting MM to PER activity
                    paramsPER, paramVarPER = curve_fit(MM, repDF["SubConcen"],
                                                       repDF["PER activity"],
                                                       bounds=(0, np.inf))
                    PERcol = "PER {0:d}".format(repPar)
                    VmaxDF.loc[sampleIndex, PERcol] = paramsPER[0]
                    KmDF.loc[sampleIndex, PERcol] = paramVarPER[1]
                except RuntimeError:
                    print("Can't fit {0:}, PPO replicate {1:}".format(sample,
                                                                      repPar))
                    VmaxDF.loc[sampleIndex, PERcol] = "can't fit"
                    KmDF.loc[sampleIndex, PERcol] = "can't fit"

    # Reformatting eventual hydrolase parameters dataframe to have columns ID,
    # Enzyme, Vmax, and Km. Oxidase dataframe is also reformatted like this,
    # but their value for "Enzyme" contains both the enzyme name and replicate
    # number, which will be processed in the following if statement.
    VmaxDF = pd.melt(VmaxDF, id_vars="ID", var_name="Enzyme",
                     value_name="Vmax")
    KmDF = pd.melt(KmDF, id_vars="ID", var_name="Enzyme", value_name="Km")
    paramsDF = pd.merge(left=VmaxDF, right=KmDF, how="inner",
                        on=["ID", "Enzyme"])

    # Separates values in the "Enzyme" column of the oxidase parameters
    # dataframe into enzyme and replicate, with the enzyme name going back
    # to the original "Enzyme" column and the replicate number being put in
    # a new column.
    if enzymeType == "O":
        paramsDF["Enzyme"] = paramsDF["Enzyme"].str.split(" ")
        for index, row in paramsDF.iterrows():
            enzymeList = row["Enzyme"]
            paramsDF.loc[index, "Replicate"] = enzymeList[1]
            paramsDF.loc[index, "Enzyme"] = enzymeList[0]
    paramsDF = paramsDF.sort_values(by=["ID"])
    return paramsDF


def plotRegress(data, params, enzymeType, processInstance, timepoint,
                saveFolder=None):
    """
    Plots the cleaned enzyme activity after removing substrate inhibition and
    setting negative values to 0 along with predicted enzyme activity from
    their fitted Michaelis-Menten parameters.

    Parameters
    ----------
    data : Pandas dataframe
        Dataframe containing either hydrolase activity or oxidase activity of a
        particular time point. The dataframe should not contain both types
        of enzymes, only 1 type.
    params : Pandas dataframe
        Dataframe containing Michaelis-Menten parameters of hydrolases or
        oxidases for a particular time point.
    enzymeType : string
        Values are either "H" for hydrolases or "O" for oxidases. Dataframes
        of hydrolase and oxidase activity were formatted differently, so they
        will be treated differently in this function
    processInstance : int
        The number of times that the calculated enzyme activity had been
        processed/cleaned. The 1st time that it is cleaned includes the simple
        replacement of negative activity values with 0 and removing substrate
        inhibition.
    timepoint : int
        The time point for when the litter bags were sampled and removed from
        the field. The only purpose of this parameter is to deal with sample
        47RRX during time point 5, when its oxidase activity was assayed twice.
    saveFolder : Path object, string, or None (by default)
        If it's a Path object or string, then it is the folder to which the
        figures will be saved. By default it is None, and so the figures would
        not be saved.

    Returns
    -------
    None.

    """
    # Removing sample 47RRX from oxidase T5 as it was assayed twice, to process
    # it separately. The dates for these 2 samples are 190125 & 190222
    if enzymeType == "O" and timepoint == 5:
        samp47 = data[data["ID"] == "47RRX"]
        data = data[data["ID"] != "47RRX"]
    # Obtaining hydrolase enzyme names if data is of hydrolase
    elif enzymeType == "H":
        hydroEnzymes = data.groupby("Enzyme")["Enzyme"].count.index.tolist()
    samples = data.groupby("ID")["ID"].count.index.tolist()

    # Just in case the saveFolder is provided as an argument but isn't a Path
    # object
    if saveFolder is not None and type(saveFolder) == str:
        saveFolder = Path(saveFolder)

    # Plotting actual calculated enzyme activity and predicted activity based
    # on fitted Michaelis-Menten parameters for hydrolase and oxidase (that
    # doesn't include 47RRX if it's T5)
    for sample in samples:
        sampleData = data[data["ID"] == sample]
        sampleParams = params[params["ID"] == sample]
        date = sampleData.groupby("Assay date")["Assay date"].count.index
        date = date.tolist()[0]
        figTitle = "{0}, T{1}, {2}, processed {3}".format(sample, timepoint,
                                                          date,
                                                          processInstance)
        py.figure(figTitle, (25, 15))

        # If data and parameters provided are for hydrolase, then this plots
        # predicted & actual activities of hydrolase
        if enzymeType == "H":
            for i in range(len(hydroEnzymes)):
                enzyme = hydroEnzymes[i]
                paramValues = sampleParams[sampleParams["Enzyme"] == enzyme]
                VmaxHydro = paramValues["Vmax"]
                KmHydro = paramValues["Km"]
                if VmaxHydro != "can't fit":
                    # Creating subplot for this enzyme
                    py.subplot(2, 3, i + 1)
                    py.title(enzyme)
                    py.xlabel("Substrate concentration (micromolar)")
                    py.ylabel("Normalized enzyme activity (micromole/[g*hr])")

                    # Plotting this enzyme's calculated (actual) activity
                    enzymeDF = sampleData[sampleData["Enzyme"] == enzyme]
                    py.plot("SubConcen", "Activity", "o-", data=enzymeDF,
                            label="Actual values")

                    # Estimating and plotting activity based on parameters
                    maxConcenInd = enzymeDF["SubConcen"].idxmax()
                    maxConcen = enzymeDF.loc[maxConcenInd, "SubConcen"]
                    minConcenInd = enzymeDF["SubConcen"].idxmin()
                    minConcen = enzymeDF.loc[minConcenInd, "SubConcen"]
                    Sregress = np.linspace(minConcen, maxConcen)
                    Vhat = MM(Sregress, VmaxHydro, KmHydro)
                    py.plot(Sregress, Vhat, "-r", label="Michaelis-Menten")
                    py.legend()

        # If the data and parameters provided are for oxidase, then this plots
        # predicted and actual activities of oxidases (that doesn't include
        # T5 47RRX)
        elif enzymeType == "O":
            for i in range(2):
                # Obtaining calculated enzyme activity and parameters from this
                # replicate
                replicate = i + 1
                repData = sampleData[sampleData["Replicate"] == replicate]
                repParam = sampleParams[sampleParams["Replicate"] == replicate]

                # Generating substrate concentrations that will be used to
                # estimate enzyme activity based on MM parameters
                maxConcenInd = repData["SubConcen"].idxmax()
                maxConcen = repData.loc[maxConcenInd, "SubConcen"]
                minConcenInd = repData["SubConcen"].idxmin()
                minConcen = repData.loc[minConcenInd, "SubConcen"]
                Sregress = np.linspace(minConcen, maxConcen)

                # Obtaining PPO parameters to check if PPO activity can be
                # plotted
                PPOparam = repParam[repParam["Enzyme"] == "PPO"]
                VmaxPPO = PPOparam["Vmax"]
                KmPPO = PPOparam["Km"]
                if VmaxPPO != "can't fit":
                    # Creating PPO subplot
                    py.subplot(2, 2, (i*2) + 1)
                    subplotTitle = "PPO replicate {0}".format(replicate)
                    py.title(subplotTitle)
                    py.xlabel("Substrate concentration (micromolar)")
                    py.ylabel("Normalized enzyme activity (micromole/[g*hr])")

                    # Plotting PPO actual activities for this replicate
                    py.plot("SubConcen", "PPO activity", "o-", data=repData,
                            label="Actual values")

                    # Estimating & plotting PPO activities based on parameters
                    # for this replicate
                    VhatPPO = MM(Sregress, VmaxPPO, KmPPO)
                    py.plot(Sregress, VhatPPO, "-r", label="Michaelis-Menten")
                    py.legend()

                # Obtaining PER parameters to check if PER activity can be
                # plotted
                PERparam = repParam[repParam["Enzyme"] == "PER"]
                VmaxPER = PERparam["Vmax"]
                KmPER = PERparam["Km"]
                if VmaxPER != "can't fit":
                    # Creating PER subplot
                    py.subplot(2, 2, (i*2) + 2)
                    subplotTitle = "PER replicate {0}".format(replicate)
                    py.title(subplotTitle)
                    py.xlabel("Substrate concentration (micromolar)")
                    py.ylabel("Normalized enzyme activity (micromole/[g*hr])")

                    # Plotting actual PER activity
                    py.plot("SubConcen", "PER activity", "o-", data=repData,
                            label="Actual values")

                    # Estimating and plotting PER activity based on parameters
                    VhatPER = MM(Sregress, VmaxPER, KmPER)
                    py.plot(Sregress, VhatPER, "-r", label="Michaelis-Menten")
                    py.legend()

        # Saving figure if the user intentionally provided a folder to save
        if saveFolder is not None:
            savePath = saveFolder/figTitle/".png"
            py.savefig(savePath)

    # If data is T5 oxidase, then this plots sample 47RRX
    if enzymeType == "O" and timepoint == 5:
        sampleParams = params[params["ID"] == "47RRX"]
        dates = samp47.groupby("Assay date")["Assay date"].count.index.tolist()
        for date in dates:
            dateData = samp47[samp47["Assay date"] == date]
            figTitle = "{0}, T{2}, {1}, processed {3}".format(sample, date,
                                                              timepoint,
                                                              processInstance)

            # Obtaining parameters for a specific date. For a specific date,
            # there is a set of parameter each for 2 enzymes and 2 replicates,
            # so for a specific date, there are 4 sets of parameters in total.
            if date == "190125":
                dateParams = sampleParams[sampleParams["Replicates"] <= 2]
            elif date == "190222":
                dateParams = sampleParams[sampleParams["Replicates"] >= 3]
            py.figure(figTitle, (25, 15))

            # Obtaining replicates from dateParams
            replicates = dateParams.groupby("Replicates")["Replicates"].count
            replicates = replicates.index.tolist()
            for index, replicate in replicates.iteritems():
                # Since in the T5 oxidase calculated enzyme activities data
                # frame, the replicates are still put as 1 and 2 for sample
                # 47RRX when, respectively, they are 3 and 4 in the parameters
                # dataframe, the following if statement converts the replicate
                # numbers from dateParams to replicate numbers used in the
                # calculated enzyme activities dataframe
                if replicate >= 3:
                    dataRep = replicate - 2  # replicate used in the dataframe
                    # of calculated enzyme activities

                    repData = dateData[dateData["Replicate"] == dataRep]
                    # the dataframe of calculated activities that contains this
                    # replicate

                # Generating substrate concentrations that will be used to
                # estimate enzyme activity based on MM parameters
                maxConcenInd = repData["SubConcen"].idxmax()
                maxConcen = repData.loc[maxConcenInd, "SubConcen"]
                minConcenInd = repData["SubConcen"].idxmin()
                minConcen = repData.loc[minConcenInd, "SubConcen"]
                Sregress = np.linspace(minConcen, maxConcen)

                # Obtaining parameters for this replicate from this date
                repParam = dateParams[dateParams["Replicate"] == replicate]

                # Obtaining PPO parameters for this replicate
                PPOparam = repParam[repParam["Enzyme"] == "PPO"]
                VmaxPPO = PPOparam["Vmax"]
                KmPPO = PPOparam["Km"]
                if VmaxPPO != "can't fit":
                    # Creating subplot for PPO
                    py.subplot(2, 2, replicate*2 - 1)
                    subplotTitle = "PPO replicate {0}".format(replicate)
                    py.title(subplotTitle)
                    py.xlabel("Substrate concentration (micromolar)")
                    py.ylabel("Normalized enzyme activity (micromole/[g*hr])")

                    # Plotting calculated PPO activity at this replicate
                    py.plot("SubConcen", "PPO activity", "o-", data=repData,
                            label="Actual values")

                    # Estimating and plotting PPO activity based on parameters
                    VhatPPO = MM(Sregress, VmaxPPO, KmPPO)
                    py.plot(Sregress, VhatPPO, "-r", label="Michaelis-Menten")
                    py.legend()

                # Obtaining PER parameters for this replicate to check if PER
                # activity can be estimated and actual & estimated activity
                # can be plotted.
                PERparam = repParam[repParam["Enzyme"] == "PER"]
                VmaxPER = PERparam["Vmax"]
                KmPER = PERparam["Km"]
                if VmaxPER != "can't fit":
                    # Creating PER subplot
                    py.subplot(2, 2, replicate*2)
                    subplotTitle = "PER replicate {0}".format(replicate)
                    py.title(subplotTitle)
                    py.xlabel("Substrate concentration (micromolar)")
                    py.ylabel("Normalized enzyme activity (micromole/[g*hr])")

                    # Plotting calculated PER activity
                    py.plot("SubConcen", "PER activity", "o-", data=repData,
                            label="Actual values")

                    # Estimating & plotting PER activity using MM parameters
                    VhatPER = MM(Sregress, VmaxPER, KmPER)
                    py.plot(Sregress, VhatPER, "-r", label="Michaelis-Menten")
                    py.legend()
            if saveFolder is not None:
                savePath = saveFolder/figTitle/".png"
                py.savefig(savePath)
    return
