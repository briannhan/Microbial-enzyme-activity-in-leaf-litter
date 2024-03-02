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
Vmax = (pd.read_excel(VmaxPath, dtype={"timePoint": "string"})
        .rename(columns={"Enzyme": "dependent", "Vmax": "data"})
        )
VmaxPartialEta2cols = {"Enzyme": "dependent", "Precipitation": "Precip",
                       "Vegetation x Precipitation": "interaction"}
VmaxPartialEta2 = (pd.read_excel(partialEta2path, "Vmax",
                                 usecols=[0, 1, 2, 3, 6])
                   .rename(columns=VmaxPartialEta2cols)
                   )
litterChemPartialEta2cols = {"functionalGroup": "dependent",
                             "Precipitation": "Precip",
                             "Vegetation x Precipitation": "interaction"}
litterChemPartialEta2 = (pd.read_excel(partialEta2path, "litterChemistry",
                                       usecols=[0, 1, 2, 3, 6])
                         .rename(columns=litterChemPartialEta2cols)
                         .query("dependent != ['carboEster', 'amide']")
                         )
VmaxLme = pd.read_excel(lmePath, "Vmax")
litterChemLme = pd.read_excel(lmePath, "litterChem")
litterChemLme.loc[litterChemLme.transformation.isna(), "transformation"] = None
independentVars = VmaxLme.columns.tolist()[2:]
mergeCols = VmaxLme.columns.tolist()[2:]
mergeCols.append("dependent")

# Filtering out the insignificant dependent variables from the partial eta^2
# and linear mixed effect model dataframes
VmaxPartialEta2 = VmaxPartialEta2.loc[VmaxPartialEta2.Vegetation.notna() | VmaxPartialEta2.Precip.notna() | VmaxPartialEta2.interaction.notna()]
litterChemPartialEta2 = litterChemPartialEta2.loc[litterChemPartialEta2.Vegetation.notna() | litterChemPartialEta2.Precip.notna() | litterChemPartialEta2.interaction.notna()]
VmaxLme = VmaxLme.loc[VmaxLme.Vegetation.notna() | VmaxLme.Precip.notna() | VmaxLme.interaction.notna()]
litterChemLme = litterChemLme.loc[litterChemLme.Vegetation.notna() | litterChemLme.Precip.notna() | litterChemLme.interaction.notna()]
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


def Tukey(data, partialEta2, lme, defaultTransformation=None):
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
    defaultTransformation : str, optional
        The original transformation put onto the data for the initial
        ANOVA + Tukey's analysis. The default is None.

    Returns
    -------
    None. Exports Tukey's results if they are performed.

    """
    partialEta2copy = partialEta2.copy()
    lmeCopy = lme.copy()
    for column in independentVars:
        partialEta2copy.loc[partialEta2copy[column].notna(), column] = "*"
        lmeCopy.loc[lmeCopy[column].notna(), column] = "*"
    merge = lmeCopy.merge(partialEta2copy, "outer", mergeCols, indicator=True)
    for index, row in merge.iterrows():
        lmeTrans = row["transformation"]
        dependent = row["dependent"]
        subsetData = data.query("dependent == @dependent")
        if lmeTrans != defaultTransformation:
            for var in independentVars:
                if var == "interaction":
                    subsetData["interaction"] = subsetData["Vegetation"] + " x " + subsetData["Precip"]
                posthoc = pairwise_tukeyhsd(subsetData.data, subsetData[var])
                posthoc = posthoc.summary().data
                
    return
