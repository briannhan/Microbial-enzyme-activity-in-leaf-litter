# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 20:35:38 2024

@author: Brian Chung
This script performs Tukey's post-hoc on linear mixed effect model results. For
now, linear mixed models have only been ran on Vmax and FTIR bands. While I do
have a CAZyme gene dataset, Ashish said that he will give me a better dataset
to work on that, I presume, will look at %reads instead of %CAZyme gene
proportions.
"""
import os
from pathlib import Path
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from compactletterdisplay.cld_calculator import compact_letter_display as cld

# Defining paths to the relevant directories
repository = Path(os.getcwd())
statsFolder = repository/"Statistical analyses"
cldFolder = statsFolder/"Compact letter displays"
enzymeCLD = cldFolder/"Enzyme Vmax"
litterChemCLD = cldFolder/"Litter chemistry"
TukeyPostHocFolder = statsFolder/"Tukey posthoc"
enzymeFolder = repository/"Enzyme activity data"
litterChemFolder = repository/"Litter chemistry"

"""Importing the relevant files

Litter chemistry FTIR.xlsx
- Data, 1st tab
- no need to separate out the concatenated timepoint & plot ID

Vmax.xlsx

Partial eta-squared.xlsx
- Vmax, 1st tab
- litterChemistry, 2nd tab

Linear mixed models, Cohen's D for main effects.xlsx
- Vmax, 2nd tab
- litterChem, 3rd tab

Let's rename the columns across all these dataframes so that they are uniform.
The results dataframes (Partial eta-squared & linear mixed models) will have
their independent variable columns rename to match that of the FTIR and Vmax
data. The 'Enzyme' columns in the Vmax data and its partial eta-squared results
and the 'functionalGroup' columns in the FTIR data and its partial eta-squared
results will be renamed to 'dependent'. From the partial eta-squared results,
the only interaction that will be kept will be 'Vegetation x Precipitation'.
This column will be renamed to 'interaction'
"""
litterChemPath = litterChemFolder/"Litter chemistry FTIR.xlsx"
VmaxPath = enzymeFolder/"Vmax.xlsx"
partialEta2path = statsFolder/"Partial eta-squared.xlsx"
lmePath = statsFolder/"Linear mixed models, Cohen's D for main effects.xlsx"
litterChem = (pd.read_excel(litterChemPath, "Data")
              .rename(columns={"functionalGroup": "dependent",
                               "spectralArea": "data"})
              )
litterChem["parameter"] = litterChem["dependent"]
Vmax = (pd.read_excel(VmaxPath, dtype={"timePoint": "string"})
        .rename(columns={"Enzyme": "dependent", "Vmax": "data"})
        )
Vmax["parameter"] = "Vmax"
VmaxPartialEta2cols = {"Enzyme": "dependent", "Precipitation": "Precip",
                       "Vegetation x Precipitation": "interaction"}
VmaxPartialEta2 = (pd.read_excel(partialEta2path, "Vmax",
                                 usecols=[0, 1, 2, 3, 6])
                   .rename(columns=VmaxPartialEta2cols)
                   )
VmaxPartialEta2["transformation"] = "log10"
VmaxPartialEta2["parameter"] = "Vmax"
litterChemPartialEta2cols = {"functionalGroup": "dependent",
                             "Precipitation": "Precip",
                             "Vegetation x Precipitation": "interaction"}
litterChemPartialEta2 = (pd.read_excel(partialEta2path, "litterChemistry",
                                       usecols=[0, 1, 2, 3, 6])
                         .rename(columns=litterChemPartialEta2cols)
                         .query("dependent != ['carboEster', 'amide']")
                         )
litterChemPartialEta2["transformation"] = None
litterChemPartialEta2["parameter"] = litterChemPartialEta2["dependent"]
VmaxLme = pd.read_excel(lmePath, "Vmax")
VmaxLme["parameter"] = "Vmax"
litterChemLme = pd.read_excel(lmePath, "litterChem")
litterChemLme.loc[litterChemLme.transformation.isna(), "transformation"] = None
litterChemLme["parameter"] = litterChemLme["dependent"]
independentVars = VmaxLme.columns.tolist()[2:-1]
mergeCols = VmaxLme.columns.tolist()[1:]
mergeCols.append("dependent")

# Filtering out the insignificant dependent variables from the partial eta^2
# and linear mixed effect model dataframes
# VmaxPartialEta2 = VmaxPartialEta2.loc[VmaxPartialEta2.Vegetation.notna() | VmaxPartialEta2.Precip.notna() | VmaxPartialEta2.interaction.notna()]
# litterChemPartialEta2 = litterChemPartialEta2.loc[litterChemPartialEta2.Vegetation.notna() | litterChemPartialEta2.Precip.notna() | litterChemPartialEta2.interaction.notna()]
# VmaxLme = VmaxLme.loc[VmaxLme.Vegetation.notna() | VmaxLme.Precip.notna() | VmaxLme.interaction.notna()]
# litterChemLme = litterChemLme.loc[litterChemLme.Vegetation.notna() | litterChemLme.Precip.notna() | litterChemLme.interaction.notna()]
# %%
"""Performing Tukey's post-hoc.

I will only perform Tukey's post-hoc under the following criteria
- if the transformation for Vmax under linear mixed effect models aren't log10
and if the litter chemistry data is transformed for linear mixed effect models
- Should the transformation between linear mixed effect models and the previous
analysis of ANOVA + Tukey's be the same and linear mixed effect models found a
significant effect where ANOVA + Tukey's did not

Following the 2 criteria above, I will not perform Tukey's post-hoc under the
following scenario
- if the transformation is the same between linear mixed effect models and
ANOVA + Tukey's and both linear mixed effect models and ANOVA + Tukey's found
a significant effect
"""

# Writing a function to perform the Tukey's


def Tukey(data, partialEta2, lme):
    """
    Performs Tukey's post-hoc following linear mixed effect model results using
    the conditions above. This function performs Tukey's post-hoc on all
    dependent variables in a dataset (e.g. all enzymes in the Vmax dataset, all
    FTIR bands in the litter chemistry dataset, etc) if the transformation for
    linear mixed effect models are the same as ANOVA + Tukey's and linear mixed
    effect models found a significant effect where ANOVA + Tukey's did not.

    Parameters
    ----------
    data : Pandas dataframe
        The dataset which Tukey's post-hoc will be performed upon.
    partialEta2 : Pandas dataframe
        Partial eta-squared from a previous ANOVA + Tukey's analysis of the
        dataset.
    lme : Pandas dataframe
        Linear mixed effect model results, with Cohen's D for significant
        main effects, of the dataset.

    Returns
    -------
    A dataframe of dependent variables and independent variables that have
    underwent Tukey's post-hoc.

    Exports Tukey's results if they are performed.

    """
    partialEta2copy = partialEta2.copy()
    lmeCopy = lme.copy()
    for column in independentVars:
        partialEta2copy.loc[partialEta2copy[column].notna(), column] = "*"
        partialEta2copy.loc[partialEta2copy[column].isna(), column] = "x"
        lmeCopy.loc[lmeCopy[column].notna(), column] = "*"
        lmeCopy.loc[lmeCopy[column].isna(), column] = "x"
    merge = (lmeCopy.merge(partialEta2copy, "left", mergeCols, indicator=True)
             .query("_merge == 'left_only'")
             )
    merge = merge.loc[merge.Vegetation.isin(["*"]) | merge.Precip.isin(["*"]) | merge.interaction.isin(["*"])]

    # For each enzyme, functional group, or whatever dependent variable the new
    # metagenomic dataset would be, label the independent variable(s) that
    # would get Tukey'ed
    for index, row in merge.iterrows():
        lmeTrans = row["transformation"]
        dependent = row["dependent"]
        partialEta2row = partialEta2copy.loc[partialEta2.dependent == dependent]
        partialEta2trans = partialEta2row["transformation"].tolist()[0]
        if lmeTrans != partialEta2trans:
            for var in independentVars:
                if row[var] == "*":
                    merge.loc[index, var] = "yes"
        elif lmeTrans == partialEta2trans:
            for var in independentVars:
                partialEta2rowVar = partialEta2row[var].tolist()[0]
                if row[var] == "*" and partialEta2rowVar == "x":
                    merge.loc[index, var] = "yes"

    # Actually performing Tukey's post-hoc
    merge = merge.loc[merge.Vegetation.isin(["yes"]) | merge.Precip.isin(["yes"]) | merge.interaction.isin(["yes"])]
    merge = merge.drop(columns=["timePoint", "_merge"])
    for column in independentVars:
        merge.loc[merge[column] != "yes", column] = None
    for index, row in merge.iterrows():
        dependent = row["dependent"]
        subsetData = data.query("dependent == @dependent")
        subsetData["interaction"] = subsetData["Vegetation"] + " x " + subsetData["Precip"]
        for var in independentVars:
            if row[var] == "yes":
                posthoc = pairwise_tukeyhsd(subsetData.data, subsetData[var])
                posthoc = posthoc.summary().data
                parameter = subsetData.parameter.tolist()[0]
                transformation = row["transformation"]
                resultsName = "{0}, {1}, {2}.xlsx".format(parameter, var,
                                                          transformation)
                specificPostHocFolder = TukeyPostHocFolder/dependent
                resultsPath = specificPostHocFolder/resultsName
                TukeyDF = pd.DataFrame(posthoc[1:], columns=posthoc[0])
                if os.path.exists(resultsPath) is False:
                    TukeyDF.to_excel(resultsPath, index=False)
    return merge


def significantTukey(tukeyTestedDF, cldFolder):
    """
    This function looks through the dependent - independent variable
    combinations in a dataset that have undergone Tukey post-hoc after running
    linear mixed effect models and see which combinations are significantly
    different under Tukey's post-hoc.

    Parameters
    ----------
    tukeyTestedDF : Pandas dataframe
        A dataframe where each row is a dependent variable and the columns
        describe the independent variables that have undergone Tukey testing.
    cldFolder : Path object
        The path to the folder containing the compact letter displays for a
        specific dataset

    Returns
    -------
    A dataframe showing the dependent - independent variable combinations that
    are significant under Tukey's post-hoc.

    """
    # Creating compact letter displays and determining which
    # dependent - independent variable combinations are significant under
    # Tukey's post-hoc
    for index, row in tukeyTestedDF.iterrows():
        dependent = row["dependent"]
        transformation = row["transformation"]
        parameter = row["parameter"]
        specificPostHocFolder = TukeyPostHocFolder/dependent
        for var in independentVars:
            if row[var] == "yes":
                tukeyResultsName = "{0}, {1}, {2}.xlsx".format(parameter, var,
                                                               transformation)
                tukeyResultsPath = specificPostHocFolder/tukeyResultsName
                tukeyResultsDF = (pd.read_excel(tukeyResultsPath,
                                                dtype={"reject": "bool"})
                                  .query("reject == True")
                                  )
                if tukeyResultsDF.empty is False:
                    group1 = tukeyResultsDF.group1.tolist()
                    group2 = tukeyResultsDF.group2.tolist()
                    pairs = list(zip(group1, group2))
                    allGroups = list(set(group1 + group2))
                    display = cld(pairs, allGroups)
                    displayDF = pd.DataFrame({"groups": allGroups,
                                              "labels": display})
                    cldName = "{0}, {1}, {2}, Tukey labels.xlsx".format(dependent, var, transformation)
                    cldExportPath = cldFolder/cldName
                    tukeyTestedDF.loc[index, var] = "significant"
                    if os.path.exists(cldExportPath) is False:
                        displayDF.to_excel(cldExportPath, index=False)

    # Returning a dataframe showing signficant combinations
    significant = (tukeyTestedDF.copy()
                   .query("Vegetation == 'significant' | Precip == 'significant' | interaction == 'significant'")
                   )
    return significant


def updateTestResults(lme, tukeyTested, significant):
    """
    Updates the dataframe of linear mixed effect model results + Cohen's D by
    removing Cohen's D values for independent variables that are insignificant
    under Tukey's post-hoc

    Parameters
    ----------
    lme : Pandas dataframe
        A dataframe containing linear mixed effect model results and Cohen's D
        for significant independent variables (minus the 2-way interaction).
    tukeyTested : Pandas dataframe
        A dataframe containing the dependent - independent variable
        combinations that have undergone Tukey's post-hoc after linear mixed
        effect models.
    significant : Pandas dataframe
        A dataframe containing significant dependent - independent variable
        combinations after Tukey's post-hoc.

    Returns
    -------
    A dataframe of the updated dataframe of significant linear mixed effect
    model results, updated with Tukey's post-hoc.

    """
    updatedLme = lme.copy()
    for index, row in VmaxTukeyTested.iterrows():
        dependent = row["dependent"]
        for var in independentVars:
            tukeyResult = significant.loc[significant.dependent == dependent, var]
            if row[var] == "yes" and (tukeyResult.empty or tukeyResult.tolist()[0] is None):
                updatedLme.loc[updatedLme.dependent == dependent, var] = None

    for var in independentVars:
        updatedLme.loc[updatedLme[var].isna(), var] = None
    return updatedLme


# %%
VmaxTukeyTested = Tukey(Vmax, VmaxPartialEta2, VmaxLme)
VmaxTukeySignificant = significantTukey(VmaxTukeyTested, enzymeCLD)
VmaxFinal = updateTestResults(VmaxLme, VmaxTukeyTested, VmaxTukeySignificant)
