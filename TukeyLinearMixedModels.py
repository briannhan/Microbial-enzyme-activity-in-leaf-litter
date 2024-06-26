# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 20:35:38 2024

@author: Brian Chung
This script performs Tukey's post-hoc on linear mixed effect model results. For
now, linear mixed models have only been ran on Vmax and FTIR bands. While I do
have a CAZyme gene dataset, Ashish said that he will give me a better dataset
to work on that, I presume, will look at %reads instead of %CAZyme gene
proportions.

And now he had given me that dataset. It sums up the number of CAZyme genes for
a putative substrate per million reads. I'm assuming that this is essentially
a weighted sum of the number of unique genes multiplied by their copies, so
this measure is essentially the number of gene copies for a putative substrate
per million reads. So an increase in this number could either mean that there
are an increasing number of unique genes (e.g. an increase of different taxa
or increasing evolutionary rates) or an increasing abundance of a specific
taxon.
"""
import os
from pathlib import Path
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from compactletterdisplay.cld_calculator import compact_letter_display as cld
import numpy as np

# Defining paths to the relevant directories
repository = Path(os.getcwd())
statsFolder = repository/"Statistical analyses"
cldFolder = statsFolder/"Compact letter displays"
enzymeCLD = cldFolder/"Enzyme Vmax"
litterChemCLD = cldFolder/"Litter chemistry"
metagenomicCLD = cldFolder/"CAZyme gene per million reads"
TukeyPostHocFolder = statsFolder/"Tukey posthoc"
enzymeFolder = repository/"Enzyme activity data"
litterChemFolder = repository/"Litter chemistry"
VmaxDataPath = repository/"Enzyme activity data"/"Vmax.xlsx"
litterChemDataPath = repository/"Litter chemistry"/"Litter chemistry FTIR.xlsx"
metagenomicFolder = repository/"CAZyme metagenomic data"
metagenomicDataPath = metagenomicFolder/"CAZyme gene counts.csv"

"""I'm going to import the following files

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
- CAZyme metagenomic, 4th tab

CAZyme gene counts.csv

Let's rename the columns across all these dataframes so that they are uniform.
The results dataframes (Partial eta-squared & linear mixed models) will have
their independent variable columns rename to match that of the FTIR and Vmax
data. The 'Enzyme' columns in the Vmax data and its partial eta-squared results
and the 'functionalGroup' columns in the FTIR data and its partial eta-squared
results will be renamed to 'dependent'. From the partial eta-squared results,
the only interaction that will be kept will be 'Vegetation x Precipitation'.
This column will be renamed to 'interaction'
"""
lmePath = statsFolder/"Linear mixed models, Cohen's D for main effects.xlsx"
partialEta2path = statsFolder/"Partial eta-squared.xlsx"

"""Importing litter chemistry data, factorial ANOVA + Tukey's results for
litter chemistry, and linear mixed effect model results for litter chemistry"""
litterChemPath = litterChemFolder/"Litter chemistry FTIR.xlsx"
litterChem = (pd.read_excel(litterChemPath, "Data")
              .rename(columns={"functionalGroup": "dependent",
                               "spectralArea": "data"})
              )
litterChem["parameter"] = litterChem["dependent"]
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
litterChemLme = pd.read_excel(lmePath, "litterChem")
litterChemLme.loc[litterChemLme.transformation.isna(), "transformation"] = None
litterChemLme["parameter"] = litterChemLme["dependent"]

"""Importing Vmax data, factorial ANOVA + Tukey's for Vmax, and linear mixed
effect model results for Vmax"""
VmaxPath = enzymeFolder/"Vmax.xlsx"
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
VmaxLme = pd.read_excel(lmePath, "Vmax")
VmaxLme["parameter"] = "Vmax"
independentVars = VmaxLme.columns.tolist()[2:-1]
mergeCols = VmaxLme.columns.tolist()[1:]
mergeCols.append("dependent")

"""Importing the new metagenomic dataset and linear mixed effect model
results for the dataset"""
metagenomic = (pd.read_csv(metagenomicDataPath)
               .rename(columns={"substrate": "dependent",
                                "genesPerMillionReads": "data"})
               )
metagenomic["dependent"] = metagenomic.dependent.str.capitalize()
metagenomic["parameter"] = "Genes per million reads"
metagenomicLme = pd.read_excel(lmePath, "CAZyme metagenomic")
metagenomicLme["parameter"] = "Genes per million reads"
metagenomicLme.dependent = metagenomicLme.dependent.str.capitalize()
metagenomicLme.loc[metagenomicLme.transformation.isna(), "transformation"] = None

"""Creating a fake dataframe of partial-eta2 for the metagenomic data to
facilitate running the functions below"""
fakeNull = metagenomicLme.shape[0]*[None]
fakeResults = {"dependent": metagenomicLme.dependent, "timePoint": fakeNull,
               "Vegetation": fakeNull, "Precip": fakeNull,
               "interaction": fakeNull, "transformation": fakeNull,
               "parameter": metagenomicLme.parameter}
metagenomicPartialEta2 = pd.DataFrame(fakeResults)
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
    # Isolating the dependent - independent variable combinations that have
    # not undergone Tukey's from the previous ANOVA + Tukey's analysis, so
    # that there are less combinations to undergo Tukey's post-hoc this time
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


def generateCLD(tukeyResults):
    """
    A custom function to generate compact letter displays from Tukey's results.
    The imported function only looks at significant pairs, ignoring groups that
    are not significantly different from any groups. For each group that is
    not significantly different from any other groups, this function will take
    the compact letter displays generated by the imported function for the
    other groups and then assign their letters to this group. This function
    generates a dataframe of compact letter displays for all available groups.

    Parameters
    ----------
    tukeyResults : Pandas dataframe
        Tukey's Honest Significant Difference test results for a specific
        dependent - independent variable combination.

    Returns
    -------
    A dataframe representing the compact letter display for all groups.

    """
    sigTukeyDF = tukeyResults.query("reject == True")
    if sigTukeyDF.empty is False:
        sigGroup1 = sigTukeyDF.group1.tolist()
        sigGroup2 = sigTukeyDF.group2.tolist()
        sigPairs = list(zip(sigGroup1, sigGroup2))
        sigGroups = list(set(sigGroup1 + sigGroup2))
        display = cld(sigPairs, sigGroups)
        displayDF = pd.DataFrame({"groups": sigGroups, "labels": display})
        group1 = tukeyResults.group1.tolist()
        group2 = tukeyResults.group2.tolist()
        allGroups = set(group1 + group2)
        missingGroups = list(allGroups - set(sigGroups))

        # This if statement is to handle groups that are not significantly
        # different from any groups
        if len(missingGroups) > 0:
            additionalCLD = []
            for group in missingGroups:
                groupResults = tukeyResults.query("(group1 == @group or group2 == @group) and reject == False")
                similarGroup1 = set(groupResults.group1.tolist())
                similarGroup2 = set(groupResults.group2.tolist())
                similarGroups = similarGroup1 | similarGroup2 - {group}
                similarGroupCLD = displayDF.loc[displayDF.groups.isin(similarGroups), "labels"]
                similarGroupCLD = similarGroupCLD.drop_duplicates().tolist()
                similarGroupCLD.sort()
                similarGroupCLD = "".join(similarGroupCLD)
                additionalCLD.append([group, similarGroupCLD])
            additionalCLD = pd.DataFrame(additionalCLD, columns=["groups", "labels"])
            displayDF = pd.concat([displayDF, additionalCLD])
        return displayDF
    elif sigTukeyDF.empty:
        return


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
    significant = tukeyTestedDF.copy()
    for index, row in significant.iterrows():
        dependent = row["dependent"]
        transformation = row["transformation"]
        parameter = row["parameter"]
        specificPostHocFolder = TukeyPostHocFolder/dependent
        for var in independentVars:
            if row[var] == "yes":
                tukeyResultsName = "{0}, {1}, {2}.xlsx".format(parameter, var, transformation)
                tukeyResultsPath = specificPostHocFolder/tukeyResultsName
                tukeyResultsDF = pd.read_excel(tukeyResultsPath,
                                               dtype={"reject": "bool"})
                displayDF = generateCLD(tukeyResultsDF)

                if displayDF is not None:
                    # Exporting the compact letter display
                    cldName = "{0}, {1}, {2}, Tukey labels.xlsx".format(dependent, var, transformation)
                    cldExportPath = cldFolder/cldName
                    significant.loc[index, var] = "significant"
                    if os.path.exists(cldExportPath) is False:
                        displayDF.to_excel(cldExportPath, index=False)

    # Returning a dataframe showing signficant combinations
    for var in independentVars:
        significant.loc[significant[var] != "significant", var] = None
    significant = significant.query("Vegetation == 'significant' | Precip == 'significant' | interaction == 'significant'")

    return significant


def updateTestResults(lme, tukeyTested, significant, partialEta2, data):
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
    partialEta2 : Pandas dataframe
        A dataframe of partial eta squared values from the original analysis
        of 3-way ANOVA + Tukey's post-hoc
    data : Pandas dataframe
        A dataframe of the dataset

    Returns
    -------
    A dataframe of the updated dataframe of significant linear mixed effect
    model results, updated with Tukey's post-hoc.

    """
    updatedLme = lme.copy()

    # Updating the linear mixed effect model results with Tukey's post-hoc
    if significant.empty is False:
        '''For a dependent variable in which an independent variable had been
        found to be significant, this for loop updates the linear mixed effect
        model results to declare this combination as significant. However,
        this for loop ignores dependent variables that have been tested but
        have not been found significant.'''
        for index, row in significant.iterrows():
            dependent = row["dependent"]
            for var in independentVars:
                if row[var] == "significant":
                    if var == "interaction":
                        updatedLme.loc[updatedLme.dependent == dependent, var] = "*"
                elif row[var] is None:
                    updatedLme.loc[updatedLme.dependent == dependent, var] = None

        '''For dependent variables which had been tested under Tukey's
        post-hoc but whose independent variables were found to be
        insignificant, this updates these dependent variables to reflect their
        insignificance'''
        significantDependent = set(significant.dependent.tolist())
        testedDependent = set(tukeyTested.dependent.tolist())
        insignificantDependent = testedDependent - significantDependent
        for var in independentVars:
            updatedLme.loc[updatedLme.dependent.isin(insignificantDependent), var] = None
    elif significant.empty:
        for index, row in tukeyTested.iterrows():
            dependent = row["dependent"]
            for var in independentVars:
                if row[var] == "yes":
                    updatedLme.loc[updatedLme.dependent == dependent, var] = None

    for var in independentVars:
        updatedLme.loc[updatedLme[var].isna(), var] = pd.NA
        updatedLme.loc[updatedLme[var].isna(), var] = None

    # Updating the linear mixed effect model results with additional
    # results found significant under ANOVA + Tukey's
    partialEta2copy = partialEta2.copy()
    for var in independentVars:
        partialEta2copy.loc[partialEta2copy[var].notna(), var] = "*"

    if "Enzyme" in data.columns:
        parameter = "Vmax"
    elif "functionalGroup" in data.columns:
        parameter = "spectralArea"

    for index, row in updatedLme.iterrows():
        lmeTrans = row.transformation
        dependent = row.dependent
        subsetData = data.query("dependent == @dependent")
        if lmeTrans == "log10":
            subsetData["data"] = np.log10(subsetData["data"])
        elif lmeTrans == "reciprocal":
            subsetData["data"] = 1/subsetData["data"]
        partialEta2row = partialEta2copy.loc[partialEta2.dependent == dependent]
        partialEta2trans = partialEta2row["transformation"].tolist()[0]
        if lmeTrans == partialEta2trans:
            for var in independentVars:
                if partialEta2row[var].tolist()[0] == "*" and (isinstance(row[var], float) is False or isinstance(row[var], str) is False):
                    if var != "interaction":
                        summary = subsetData.groupby(var)["data"].agg(["mean", "var", "size"])
                        mean1 = summary.iloc[0, 0]
                        var1 = summary.iloc[0, 1]
                        size1 = summary.iloc[0, 2]
                        mean2 = summary.iloc[1, 0]
                        var2 = summary.iloc[1, 1]
                        size2 = summary.iloc[1, 2]
                        pooledVar = (((size1 - 1)*var1) + ((size2 - 1)*var2))/(size1 + size2 - 2)
                        pooledStd = pooledVar**(1/2)
                        CohenD = (mean1 - mean2)/pooledStd
                        CohenD = np.abs(CohenD)
                        updatedLme.loc[index, var] = CohenD
                    elif var == "interaction":
                        updatedLme.loc[index, var] = "*"
    return updatedLme


# %%

"""Performing Tukey's post-hoc, generating compact letter displays, and
updating test results for most of these dependent - independent variable
combinations"""
VmaxTukeyTested = Tukey(Vmax, VmaxPartialEta2, VmaxLme)
VmaxTukeySignificant = significantTukey(VmaxTukeyTested, enzymeCLD)
VmaxFinal = updateTestResults(VmaxLme, VmaxTukeyTested, VmaxTukeySignificant,
                              VmaxPartialEta2, Vmax)

litterChemTukeyTested = Tukey(litterChem, litterChemPartialEta2, litterChemLme)
litterChemTukeySignificant = significantTukey(litterChemTukeyTested,
                                              litterChemCLD)
litterChemFinal = updateTestResults(litterChemLme, litterChemTukeyTested,
                                    litterChemTukeySignificant,
                                    litterChemPartialEta2, litterChem)

metagenomicTukeyTested = Tukey(metagenomic, metagenomicPartialEta2,
                               metagenomicLme)
metagenomicSignificant = significantTukey(metagenomicTukeyTested,
                                          metagenomicCLD)
metagenomicFinal = updateTestResults(metagenomicLme, metagenomicTukeyTested,
                                     metagenomicSignificant,
                                     metagenomicPartialEta2, metagenomic)
# %%
"""The following dependent - independent variable combinations are missing
compact letter displays

From the Vmax dataset
- CBH, Vegetation (missing Tukey's)
- NAG, Vegetation (missing Tukey's)
- For these 2 combinations, I did not put them through Tukey's again because
they previously have undergone Tukey's. Their compact letter displays were
never created, likely because their time x vegetation was significant in the
sense that both enzymes consistently show a difference due to vegetation across
all time points, even with changes between the time points, so creating a
vegetation compact letter display was perhaps unnecessary

From the FTIR dataset
- carboEster2, Vegetation, Precipitation, both missing Tukey's 
- glycosidicBond, Vegetation, Precipitation, both missing Tukey's
- For these 2 dependent variables, I'm assuming that these combinations are
missing because they both have been previously found to have a veg x precip
interaction, and that I was more interested in this interaction than the
independent variables separated. This time, I judge that they're missing the
precipitation compact letter displays because I was looking at the partial
eta2 results

Performing Tukey's post-hoc on them and then creating compact letter displays
for them
"""


def manualTukey(combinations, dataset, cldFolder, finalResults):
    """
    This performs Tukey's post-hoc on some manual combinations of
    independent-dependent variables that the user can specify. For example, if
    the procedure above of performing Tukey's on linear mixed effect models
    and factorial ANOVAs had missed certain combinations that should be tested,
    then the user can run this function on those missing combinations.

    Parameters
    ----------
    combinations : List
        A list of lists. The inner lists each consist of 2 elements, the
        dependent variable and the independent variable. The outer list
        consists of these inner lists. There is no additional nesting of lists.
    dataset : Pandas dataframe
        The dataset on which the additional Tukey's will be performed.
    cldFolder : Path
        The path to the directory that holds the compact letter displays for
        a dataset.
    finalResults : Pandas dataframe
        The dataframe of linear mixed effect models + previous factorial ANOVAs
        + Tukey's analyses on a dataset. The transformations for a dependent
        variable will be obtained from this dataframe.

    Returns
    -------
    None.

    """
    for combination in combinations:
        dependent = combination[0]
        independent = combination[1]
        df = dataset.query("dependent == @dependent")
        transformation = finalResults.loc[finalResults.dependent == dependent, "transformation"].tolist()[0]
        parameter = finalResults.loc[finalResults.dependent == dependent, "parameter"].tolist()[0]
        if transformation == "log10":
            df["data"] = np.log10(df["data"])
        elif transformation == "reciprocal":
            df["data"] = 1/df["data"]
        tukeyResults = pairwise_tukeyhsd(df.data, df[independent])
        tukeyResults = tukeyResults.summary().data
        tukeyResults = pd.DataFrame(tukeyResults[1:], columns=tukeyResults[0])
        specificTukeyFolder = TukeyPostHocFolder/dependent
        tukeyFileName = "{0}, {1}, {2}.xlsx".format(parameter, independent,
                                                    transformation)
        tukeyFilePath = specificTukeyFolder/tukeyFileName
        if os.path.exists(tukeyFilePath) is False:
            tukeyResults.to_excel(tukeyFilePath, index=False)
        display = generateCLD(tukeyResults)
        if display is not None:
            cldName = "{0}, {1}, Tukey labels.xlsx".format(dependent,
                                                           independent)
            cldPath = cldFolder/cldName
            if os.path.exists(cldPath) is False:
                display.to_excel(cldPath, index=False)
    return


VmaxMissingCLD = [["CBH", "Vegetation"], ["NAG", "Vegetation"]]
litterChemMissingCLD = [["carboEster2", "Vegetation"],
                        ["carboEster2", "Precip"],
                        ["glycosidicBond", "Vegetation"],
                        ["glycosidicBond", "Precip"]]

manualTukey(VmaxMissingCLD, Vmax, enzymeCLD, VmaxFinal)
manualTukey(litterChemMissingCLD, litterChem, litterChemCLD, litterChemFinal)

"""From these additional Tukey's, we see that

carboEster2 - Vegetation, significant
carboEster2 - Precipitation, significant
glycosidicBond - Vegetation, significant
glycosidicBond - Precipitation, not significant

CBH - Vegetation, significant
NAG - Vegetation, significant

Need to update the results for litter chemistry
"""
litterChemFinal.loc[litterChemFinal.dependent == "glycosidicBond", "Precip"] = None
# %%
# Exporting the updated linear mixed effect models + Tukey's post-hoc results

# Creating a ReadMe tab for the final Excel file
note = ["This file contains linear mixed effect model results, followed by",
        "Tukey's post-hoc of significant results. For a dependent variable,",
        "if an independent variable is significant under Tukey's post-hoc",
        "then Cohen's D is calculated for this combination. If the",
        "interaction between vegetation and precipitation are significant",
        "then this interaction is labeled with an asterisk.",
        "",
        "This file includes results that were found to be significant under",
        "Tukey's post-hoc after running linear mixed effect models. It also",
        "includes results found significant under a previous ANOVA + Tukey's",
        "analysis."]
readMe = pd.DataFrame({"Note": note})

# Removing the "parameter" column from the results dataframes
VmaxFinal = VmaxFinal.drop(columns="parameter")
litterChemFinal = litterChemFinal.drop(columns="parameter")
metagenomicFinal = metagenomicFinal.drop(columns="parameter")

resultsPath = statsFolder/"Linear mixed effect models, updated with Tukey.xlsx"
if os.path.exists(resultsPath) is False:
    with pd.ExcelWriter(resultsPath) as w:
        readMe.to_excel(w, "ReadMe", index=False)
        VmaxFinal.to_excel(w, "Vmax", index=False)
        litterChemFinal.to_excel(w, "litterChem", index=False)
        metagenomicFinal.to_excel(w, "CAZyme metagenomic", index=False)
