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
value units of micromole L^-1 g^-1. I divided by 2 to take into account the
fact that half the volume of each assay well consists of the pipetted
substrate (pyrogallol) and the other half consists of the filtered
homogenate.'''
pyroConcen = (pyroHighest*subProps).tolist()
pyroConcenDF = pd.DataFrame({"PlateRow": list("ABCDEFGH"),
                             "SubConcen": pyroConcen})
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
        longSampleName = row["Long sample name"]  # splitted list
        enzymeData.loc[index, "Assay date"] = longSampleName[0]
        for part in longSampleName:
            if "X" in part:
                sampleID = part
        if "B" in longSampleName:
            enzymeData.loc[index, "ID"] = "B"
        elif "X" in sampleID:
            enzymeData.loc[index, "ID"] = sampleID
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
    # (2) Graph oxidative enzyme activity
    for sample in samples:
        sampleDF = enzymeData[enzymeData["ID"] == sample]
        vegetation = sampleDF.groupby("Vegetation")["Vegetation"].count().index
        vegetation = vegetation.tolist()
        print(vegetation)
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
