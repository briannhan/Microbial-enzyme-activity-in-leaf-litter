# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 17:48:25 2024

This script creates linear mixed models to re-analyze per Ashish's request.
"""
import pandas as pd
from pathlib import Path
import os
from scipy import stats
import numpy as np

# Importing the following class to allow for me to set certain parameters
# in linear mixed effect models, including constraining correlation to 0
from statsmodels.regression.mixed_linear_model import MixedLMParams

# Importing linear mixed effect models from formula.api since this allows
# for me to construct a model using a formula
from statsmodels.formula.api import mixedlm

# Setting paths to the relevant directories
repository = Path(os.getcwd())
enzymeFolder = repository/"Enzyme activity data"
litterChemFolder = repository/"Litter chemistry"

# Importing the data
VmaxPath = enzymeFolder/"Vmax.xlsx"
VmaxDtype = {"timePoint": "category", "ID": "category", "Vegetation": "category",
             "Precip": "category", "Enzyme": "string"}
VmaxDF = (pd.read_excel(VmaxPath, dtype=VmaxDtype)
          .drop_duplicates()
          )

litterChemPath = litterChemFolder/"Litter chemistry FTIR.xlsx"
litterChemDtype = {"id": "string", "Vegetation": "category",
                   "Precip": "category", "timePoint": "int64",
                   "functionalGroup": "category"}
litterChemDF = (pd.read_excel(litterChemPath, 0, dtype=litterChemDtype)
                .drop_duplicates()
                .astype({"timePoint": "category"})
                )
# %%
# Let's convert time from a categorical variable to a numerical variable and
# re-run mixed linear models, because the model results look ugly as shit
# with time as a categorical variable

"""Litter was sampled on August 30, 2017. Bags were deployed on September 12,
2017.

T0 sampled on November 30, 2017

T3 sampled on April 11, 2018

T5 sampled in November 2018, exact date unspecified.

T6 sampled in February 2019, exact date unspecified.

I'll convert the time to days since initial litter bag deployment on September
12, 2017. For the last 2 time points, I'll just use the 15th of the month.
"""
timeDict = {"year": [2017, 2017, 2018, 2018, 2019],
            "month": [9, 11, 4, 11, 2],
            "day": [12, 30, 11, 15, 15]}
timeDF = pd.DataFrame(timeDict)
timeSeries = pd.to_datetime(timeDF)
timeDF = (pd.DataFrame({"timePoint": ["deployment", 0, 3, 5, 6],
                        "date": timeSeries})
          .astype({"timePoint": "category"})
          )
for index, row in timeDF.iterrows():
    timeDF.loc[index, "daysSinceDeployment"] = row.date - timeDF.loc[0, "date"]
timeDF.daysSinceDeployment = timeDF.daysSinceDeployment.dt.days
VmaxDF = VmaxDF.merge(timeDF, on="timePoint")
litterChemDF = litterChemDF.merge(timeDF, on="timePoint")
# %%
# Performing linear mixed effect models for alpha glucosidase, AG

"""Let's create a model with the following parameters

Group: plot ID
Fixed effects: Vegetation, Precip, Vegetation x Precip
Random effects: timePoint, plot ID (which is inherently a random factor)
"""
AG = VmaxDF.query("Enzyme == 'AG'")
fullFormula = "Vmax ~ C(Vegetation)*C(Precip)"
randomFormula = "~0+daysSinceDeployment"
randomFormula2 = "~daysSinceDeployment"
randomFormula3 = "~daysSinceDeployment - 1"

"""Formulas 3 and 1 remove the intercept for timePoint and only estimated the
slope for timePoint. Which, might not be that good of an idea"""
# AGmodel1 = mixedlm(fullFormula, AG, randomFormula, groups="ID")
# AGmodel1results = AGmodel1.fit(method=["lbfgs"])
# AGmodel1summary = AGmodel1results.summary()

AGmodel2 = mixedlm(fullFormula, AG, randomFormula2, groups="ID")
AGmodel2results = AGmodel2.fit(method=["lbfgs"])
AGmodel2summary = AGmodel2results.summary()

# AGmodel3 = mixedlm(fullFormula, AG, randomFormula3, groups="ID")
# AGmodel3results = AGmodel3.fit(method=["lbfgs"])
# AGmodel3summary = AGmodel3results.summary()

"""Ever since changing time to a numerical variable instead of a categorical
variable, running the 1st and 3rd models using the 1st and 3rd formulas
rendered a 'Singular matrix' error. Commenting them out for now"""

"""Seems like grassland has a marginal significant effect on AG Vmax, judging
from AGmodel2results. But that's not good enough. Let's check normality of
residuals from the 2nd model, which I believe is more correct"""
AGmodel2residuals = AGmodel2results.resid
AGmodel2normality = stats.shapiro(AGmodel2residuals)

"""Model residuals highly not normal. Time to log10 transform it"""
AGlog10 = AG.copy()
AGlog10["Vmax"] = np.log10(AGlog10["Vmax"])
AGmodel4 = mixedlm(fullFormula, AGlog10, randomFormula2, groups="ID")
AGmodel4results = AGmodel4.fit(method=["lbfgs"])
AGmodel4summary = AGmodel4results.summary()
AGmodel4residuals = AGmodel4results.resid
AGmodel4normality = stats.shapiro(AGmodel4residuals)
AGmodel4randomCorrelation = -0.002/((0.635*0.196)**0.5)
"""Model residuals still not normal after log10 transformation. In addition,
there is extremely negligible correlation between the random effects.
"""

"""Natural log transforming AG Vmax and then recreating model"""
AGln = AG.copy()
AGln["Vmax"] = np.log(AGln["Vmax"])
AGmodel5 = mixedlm(fullFormula, AGln, randomFormula2, groups="ID")
AGmodel5results = AGmodel5.fit(method=["lbfgs"])
AGmodel5summary = AGmodel5results.summary()
AGmodel5residuals = AGmodel5results.resid
AGmodel5normality = stats.shapiro(AGmodel5residuals)
"""Normality is identical to if it was log10 transformed."""

"""Square root transforming and then recreating model"""
AGsquareRoot = AG.copy()
AGsquareRoot["Vmax"] = AGsquareRoot["Vmax"]**(1/2)
AGmodel6 = mixedlm(fullFormula, AGsquareRoot, randomFormula2, groups="ID")
AGmodel6results = AGmodel6.fit(method=["lbfgs"])
AGmodel6summary = AGmodel6results.summary()
AGmodel6residuals = AGmodel6results.resid
AGmodel6normality = stats.shapiro(AGmodel6residuals)
"""Normality actually worsens slightly"""

"""Reciprocal transforming AG Vmax and then recreating model"""
AGreciprocal = AG.copy()
AGreciprocal["Vmax"] = 1/AGreciprocal["Vmax"]
AGmodel7 = mixedlm(fullFormula, AGreciprocal, randomFormula2, groups="ID")
AGmodel7results = AGmodel7.fit(method=["lbfgs"])
AGmodel7summary = AGmodel7results.summary()
AGmodel7residuals = AGmodel7results.resid
AGmodel7normality = stats.shapiro(AGmodel7residuals)
"""Normality worsens by close to 1 order of magnitude by p-value"""

"""Well despite the model with the log 10-transformed AG Vmax having
uncorrelated random effects, and the correct next step would be to create a
model in which the random effects are constrained to be uncorrelated, I don't
know how to create such a model when the fixed effects are categorical
variables."""

"""In summary, log transformations are best. Looking at results from the model
created from log10 transformed Vmax, I see no effects by precipitation, no
interaction effects, and a slightly significant effect by Vegetation that I'm
going to ignore."""
resultsDFcols = ["dependent", "Vegetation", "Precip",
                 "interaction", "transformation"]
AGsignificance = ["AG", None, None, None, "log10"]
# %%
# Creating linear mixed effect models for acid phosphatase, AP

AP = VmaxDF.query("Enzyme == 'AP'")
APmodel1 = mixedlm(fullFormula, AP, randomFormula2, groups="ID")
APmodel1results = APmodel1.fit(method=["lbfgs"])
APmodel1summary = APmodel1results.summary()
"""Hmm, what do you know, there are significant effects for vegetation and
precipitation

But once I converted time to a continuous variable, the significant effects
disappeared
"""

"""Checking model residuals to see if they're likely from a normal population
"""
APmodel1residuals = APmodel1results.resid
APmodel1normality = stats.shapiro(APmodel1residuals)
"""P-value of Shapiro-Wilk is 1.5*10^-14, very small, indicating it's unlikely
that the residuals are from a normal population. Time to transform it.

New p-value with time as continuous variable = 1.97*10^-15
"""

# Log10 transforming AP and running model again
APlog10 = AP.copy()
APlog10["Vmax"] = np.log10(APlog10["Vmax"])
APmodel2 = mixedlm(fullFormula, APlog10, randomFormula2, groups="ID")
APmodel2results = APmodel2.fit(method=["lbfgs"])
APmodel2summary = APmodel2results.summary()
APmodel2residuals = APmodel2results.resid
APmodel2normality = stats.shapiro(APmodel2residuals)
"""Still not quite normal, with Shapiro-Wilk normality p-value = 0.000394

New p-value with time as continuous variable = 2.12*10^-8
"""

# Natural log transform and running model again
APln = AP.copy()
APln["Vmax"] = np.log10(APln["Vmax"])
APmodel3 = mixedlm(fullFormula, APln, randomFormula2, groups="ID")
APmodel3results = APmodel3.fit(method=["lbfgs"])
APmodel3summary = APmodel3results.summary()
APmodel3residuals = APmodel3results.resid
APmodel3normality = stats.shapiro(APmodel3residuals)
"""Normality identical to log10 transformations. I'm not gonna do natural log
transformations again in the future, then.

New p-value with time as continuous variable = 2.12*10^-8
"""

# Square root transform and running again
APsquareRoot = AP.copy()
APsquareRoot["Vmax"] = APsquareRoot["Vmax"]**(1/2)
APmodel4 = mixedlm(fullFormula, APsquareRoot, randomFormula2, groups="ID")
APmodel4results = APmodel4.fit(method=["lbfgs"])
APmodel4summary = APmodel4results.summary()
APmodel4residuals = APmodel4results.resid
APmodel4normality = stats.shapiro(APmodel4residuals)
"""Normality improves slightly but still worse than log transformations

New p-value with time as continuous variable = 1.57*10^-14
"""

# Reciprocal transform and running again
APreciprocal = AP.copy()
APreciprocal["Vmax"] = 1/APreciprocal["Vmax"]
APmodel5 = mixedlm(fullFormula, APreciprocal, randomFormula2, groups="ID")
APmodel5results = APmodel5.fit(method=["lbfgs"])
APmodel5summary = APmodel5results.summary()
APmodel5residuals = APmodel5results.resid
APmodel5normality = stats.shapiro(APmodel5residuals)
"""Normality improves, p-value = 0.0002177, lower p-value than log
transformations, so still slightly lower normality than log transformations

New p-value with time as continuous variable = 0.6767. Let's keep this
"""

# Seeing what happens if I create a model with the log10 transformed data if
# I remove the insignificant interaction
mainEffects = "Vmax ~ C(Vegetation) + C(Precip)"
APmodel6 = mixedlm(mainEffects, APlog10, randomFormula2, groups="ID")
APmodel6results = APmodel6.fit(method=["lbfgs"])
APmodel6summary = APmodel6results.summary()
APmodel6residuals = APmodel6results.resid
APmodel6normality = stats.shapiro(APmodel6residuals)
"""Nope, removing the interaction actually makes the main effects worse. They
definitely are not significant"""

"""Let's just go with a log10 transformation with the interaction in the model.
After a log10 transformation, the effects for vegetation and precipitation
disappeared.

After converting time to a continuous variable, the reciprocal transformation
is the most normal. Very slight effect with precipitation (p = 0.058).
Otherwise, results are still the same for vegetation and the interaction

Might need to do Tukey post-hoc on AP later
"""

APsignificance = ["AP", None, "-", None, "log10"]
# %%
"""Writing a function to abstract the process of creating linear mixed effect
models, transforming the underlying data, and testing the normality of the
model residuals"""


def mixedLinearModel(data, parameter, transformation=None,
                     formula=fullFormula):
    """
    Creates a linear mixed model with plot ID and time point as random effects
    and vegetation, precipitation, and their interaction as fixed effects. Then
    performs a Shapiro-Wilk test on the model residuals to check for normality.

    Parameters
    ----------
    data : Pandas dataframe
        Dataframe of the parameter of interest (e.g. Vmax of a specific enzyme,
        a specific functional group in the FTIR data, etc).
    parameter : str
        The parameter of interest to be used as the dependent variable in a
        model.
    transformation : str, optional
        The user can choose to not specify this. If the user specifies this
        parameter, then the data will be transformed prior to creating a model.
        The default is None.
    formula : str, optional
        The formula describing the fixed effects of the model. The default is
        the full formula with the interaction between vegetation and
        precipitation.

    Returns
    -------
    - the summary generated from the results instance after fitting the model
    - the results object from the Shapiro-Wilk normality test
    - the results object from the fitting of the model
    """
    if transformation == "log10":
        data = data.copy()
        data[parameter] = np.log10(data[parameter])
    elif transformation == "square root":
        data = data.copy()
        data[parameter] = data[parameter]**(1/2)
    elif transformation == "reciprocal":
        data = data.copy()
        data[parameter] = 1/data[parameter]

    formula = "{} ~ C(Vegetation)*C(Precip)".format(parameter)

    model = mixedlm(formula, data, randomFormula2, groups="ID")
    modelResults = model.fit(method=["lbfgs"])
    modelSummary = modelResults.summary()
    modelResiduals = modelResults.resid
    modelNormality = stats.shapiro(modelResiduals)
    return [modelResults, modelSummary, modelNormality]


# %%
# Creating mixed linear models for beta-glucosidase, BG
BG = VmaxDF.query("Enzyme == 'BG'")
BGmodel1 = mixedLinearModel(BG, "Vmax")
"""Shapiro-Wilk p-value after converting time = 4.527*10^-16"""

# Log 10 transforming BG and re-testing
BGmodel2 = mixedLinearModel(BG, "Vmax", "log10")
"""Shapiro-Wilk p-value after converting time = 1.3198*10^-10"""


# Square root transforming BG and re-testing
BGmodel3 = mixedLinearModel(BG, "Vmax", "square root")
"""Square root transformation has better normality than log10, although still
not quite normal

After converting time, Shapiro-Wilk p-value = 6.1116*10^-15
"""

# Reciprocal transforming BG and re-testing
BGmodel4 = mixedLinearModel(BG, "Vmax", "reciprocal")
"""Improves normality but worse than log10 and square root. Let's use results
from square root transformation

After converting time, Shapiro-Wilk p-value = 7.9298*10^-6, most normal
transformation. Let's use this
"""

"""From square root transformation, effects by vegetation and precipitation
disappeared.

From log 10 transformation, vegetation effect present, no effects from
precipitation

From reciprocal transformation after converting time, there is a significant
vegetation effect, significant precipitation effect, and no significant
interaction
"""
BGsignificance = ["BG", "***", "*", None, "reciprocal"]
# %%
# Creating mixed linear models for beta-xylosidase, BX
BX = VmaxDF.query("Enzyme == 'BX'")
BXmodel1 = mixedLinearModel(BX, "Vmax")
"""Convert time -> numerical, SW p-value = 3.829*10^-14"""

# Log 10 transforming BX and re-testing
BXmodel2 = mixedLinearModel(BX, "Vmax", "log10")
"""Improves normality, but still far from normal

Convert time -> numerical, SW p-value = 1.303*10^-05
"""

# Square root transforming BX and re-testing
BXmodel3 = mixedLinearModel(BX, "Vmax", "square root")
"""Improves normality but worse than log10

Convert time -> numerical, SW p-value = 2.0909*10^-10
"""

# Reciprocal transforming BX and re-testing
BXmodel4 = mixedLinearModel(BX, "Vmax", "reciprocal")
"""Actually worsens normality

Convert time -> numerical, SW p-value = 6.72*10^-13
"""

"""Let's just use results from log 10 transformation.
No effects by either vegetation or precipitation, no interaction either"""
BXsignificance = ["BX", None, None, None, "log10"]
# %%
# Creating mixed linear models for cellobiohydrolase, CBH
CBH = VmaxDF.query("Enzyme == 'CBH'")
CBHmodel1 = mixedLinearModel(CBH, "Vmax")
"""Shapiro-Wilk p-value = 1.11*10^-6

After converting time to numerical variable, Shapiro-Wilk p-value = 0.00020449
"""

# Log 10 transforming CBH and re-testing
CBHmodel2 = mixedLinearModel(CBH, "Vmax", "log10")
"""Shapiro-Wilk p-value = 0.6667, let's keep this transformation.
Significant vegetation effect, significant precipitation effect,
no interaction

Shapiro-Wilk p-value after converting time to continuous variable = 0.495
"""

"""From log10 transformation, significant vegetation effect, marginally
significant precipitation effect that warrants at least a Tukey. No significant
interaction."""
CBHsignificance = ["CBH", "***", "-", None, "log10"]
# %%
# Creating mixed linear models for LAP
LAP = VmaxDF.query("Enzyme == 'LAP'")
LAPmodel1 = mixedLinearModel(LAP, "Vmax")
"""Shapiro-Wilk p = 4.503*10^-12"""

# Log 10 transforming and re-testing
LAPmodel2 = mixedLinearModel(LAP, "Vmax", "log10")
"""Shapiro-Wilk p = 6.254*10^-5"""

# Square root transforming and re-testing
LAPmodel3 = mixedLinearModel(LAP, "Vmax", "square root")
"""Shapiro-Wilk p = 1.81*10^-9"""

# Reciprocal transforming and re-testing
LAPmodel4 = mixedLinearModel(LAP, "Vmax", "reciprocal")
"""Shapiro-Wilk p = 0.000354, the most normal transformation. Let's use
this"""

"""From reciprocal transformation, significant vegetation effect, no
precipitation effect, no interaction"""
LAPsignificance = ["LAP", "**", None, None, "reciprocal"]
# %%
# Creating mixed linear models for NAG
NAG = VmaxDF.query("Enzyme == 'NAG'")
NAGmodel1 = mixedLinearModel(NAG, "Vmax")
"""Shapiro-Wilk normality p = 5.189*10^-7"""

# Log 10 transforming and re-testing
NAGmodel2 = mixedLinearModel(NAG, "Vmax", "log10")
"""Shapiro-Wilk p = 0.000451"""

# Square root transforming and re-testing
NAGmodel3 = mixedLinearModel(NAG, "Vmax", "square root")
"""Shapiro-Wilk normality p = 2.778*10^-6"""

# Reciprocal transforming and re-testing
NAGmodel4 = mixedLinearModel(NAG, "Vmax", "reciprocal")
"""Normality worsens, SW p = 1.0003*10^-9"""

"""From log 10 transformation, significant Vegetation effect, no precipitation
effect, no interaction"""
NAGsignificance = ["NAG", "***", None, None, "log10"]
# %%
"""Creating mixed linear models for PPO. This is gonna be more complicated
because replicate measurements of PPO were taken for each plot and timepoint.

Fixed effects: Vegetation, Precip, their interaction
Groups: plot ID (top level group, with replicate as a nested group)
Random effects: plot ID, Time (in days since deployment)
Variance component: replicate
"""
PPO = VmaxDF.query("Enzyme == 'PPO'")
PPO.loc[PPO.Replicate == 1, "Replicate"] = "a"
PPO.loc[PPO.Replicate == 2, "Replicate"] = "b"
PPO.loc[PPO.Replicate == 3, "Replicate"] = "c"
PPO.loc[PPO.Replicate == 4, "Replicate"] = "d"
PPO.Replicate = PPO.Replicate.astype("category")
vc = {"Replicate": "C(Replicate)"}
PPOmodel1 = mixedlm(fullFormula, PPO, randomFormula2, vc, groups="ID")
PPOmodel1results = PPOmodel1.fit()
PPOmodel1summary = PPOmodel1results.summary()
PPOmodel1normality = stats.shapiro(PPOmodel1results.resid)
"""Shapiro-Wilk p = 6.512*10^-24"""

# Log 10 transforming and re-testing
PPOlog10 = PPO.copy()
PPOlog10.Vmax = np.log10(PPOlog10.Vmax)
PPOmodel2 = mixedlm(fullFormula, PPOlog10, randomFormula2, vc, groups="ID")
PPOmodel2results = PPOmodel2.fit()
PPOmodel2summary = PPOmodel2results.summary()
PPOmodel2normality = stats.shapiro(PPOmodel2results.resid)
"""Shapiro-Wilk p = 3.54*10^-13"""

# Reciprocal transforming and re-testing
PPOreciprocal = PPO.copy()
PPOreciprocal.Vmax = 1/PPOreciprocal.Vmax
PPOmodel3 = mixedlm(fullFormula, PPOreciprocal, randomFormula2, vc,
                    groups="ID")
PPOmodel3results = PPOmodel3.fit()
PPOmodel3summary = PPOmodel3results.summary()
PPOmodel3normality = stats.shapiro(PPOmodel3results.resid)
"""Shapiro-Wilk p = 4.514*10^-24"""

# Square root transforming and re-testing
PPOsquareRoot = PPO.copy()
PPOsquareRoot.Vmax = (PPOsquareRoot.Vmax)**0.5
PPOmodel4 = mixedlm(fullFormula, PPOsquareRoot, randomFormula2, vc,
                    groups="ID")
PPOmodel4results = PPOmodel4.fit()
PPOmodel4summary = PPOmodel4results.summary()
PPOmodel4normality = stats.shapiro(PPOmodel4results.resid)
"""Shapiro-Wilk p = 7.981*10^-23"""

"""From the log 10 transformation, no effect by vegetation or precipitation.
Significant 2-way interaction."""
PPOsignificance = ["PPO", None, None, "**", "log10"]

"""So now all the enzymes have been tested. Moving onto the litter chemistry
functional groups."""
# %%
"""Creating a dataframe of summary results for enzyme Vmax. This dataframe
will be used to calculate effect sizes as well, using f2 as a measure of effect
size for fixed effects (Lorah 2018), which is calculated using R2 of the full
model and the R2 of the model without a specific fixed effect. R2 will be
calculated according to Snijders and Bosker (1994), using the sum of squared
residuals of the full model and the sum of squared residuals of a null model
with only 1 intercept for the fixed effect and 1 intercept for the random
effect (which will be the plot ID). When I fitted the null model to the enzyme
Vmax data, the null model's intercept is always the mean of the dependent
variable across all plots and time points. So I'll just be calculating the
mean and using that to calculate the sum of squared residuals of the null
model.


Lorah, J. (2018). Effect size measures for multilevel models: definition,
interpretation, and TIMSS example. Large-scale Assessments in Education, 6(8).
https://doi.org/10.1186/s40536-018-0061-2

Snijders, T. A. B., & Bosker, R. J. (1994) Modeled Variance in Two-Level
Models. Sociological Methods & Research, 22(3), 342-363.
https://doi.org/10.1177/0049124194022003004
"""
enzymeResultsLists = [AGsignificance, APsignificance, BGsignificance,
                      BXsignificance, CBHsignificance, LAPsignificance,
                      NAGsignificance, PPOsignificance]
enzymeResults = (pd.DataFrame(enzymeResultsLists, columns=resultsDFcols)
                 .astype(dtype={"Vegetation": "string", "Precip": "string",
                                "interaction": "string"})
                 )
for column in enzymeResults.columns:
    enzymeResults.loc[enzymeResults[column].isna(), column] = pd.NA
independent = resultsDFcols[1:4]

"""Calculating sum of squared residuals of the full model for each enzyme,
and sum of squared residuals using the null model, which is just the mean of
the dependent variable across all plots and time points."""
enzymeModels = [AGmodel4results, APmodel2results, BGmodel4, BXmodel2,
                CBHmodel2, LAPmodel4, NAGmodel2, PPOmodel2results]
for index, row in enzymeResults.iterrows():
    transformation = row["transformation"]
    enzyme = row["dependent"]
    df = VmaxDF.query("Enzyme == @enzyme")
    if transformation == "log10":
        df["Vmax"] = np.log10(df["Vmax"])
    elif transformation == "reciprocal":
        df["Vmax"] = 1/df["Vmax"]
    fullModel = enzymeModels[index]
    if isinstance(fullModel, list):
        fullModel = fullModel[0]
    fullModelResiduals = fullModel.resid
    fullModelSSR = np.sum(fullModelResiduals**2)
    enzymeResults.loc[index, "fullModelSSR"] = fullModelSSR
    nullModelSSR = (df["Vmax"].var())*(df.shape[0] - 1)
    enzymeResults.loc[index, "nullModelSSR"] = nullModelSSR
enzymeResults["fullModelR2"] = 1 - (enzymeResults.fullModelSSR/enzymeResults.nullModelSSR)

"""Fitting models without the significant fixed effects to calculate R2 of
models without significant effects. These R2 will later be used with the R2 of
the full models to calculate f2 of the significant fixed effects"""
modelsNoSignificantFixed = []
for index, row in enzymeResults.iterrows():
    transformation = row["transformation"]
    enzyme = row["dependent"]
    df = VmaxDF.query("Enzyme == @enzyme")
    if transformation == "log10":
        df["Vmax"] = np.log10(df["Vmax"])
    elif transformation == "reciprocal":
        df["Vmax"] = 1/df["Vmax"]
    # print(row.shape)
    # Fitting the models without the significant fixed effects
    modelFits = [enzyme]
    if type(row["Vegetation"]) is str:
        noVegFormula = "Vmax ~ C(Precip) + C(Vegetation):C(Precip)"
        noVegModel = mixedlm(noVegFormula, df, randomFormula2, groups="ID")
        noVegResults = noVegModel.fit(method=["lbfgs"])
        noVegResiduals = noVegResults.resid
        enzymeResults.loc[index, "noVegSSR"] = np.sum(noVegResiduals**2)
        modelFits.append(noVegResults)
    else:
        modelFits.append(pd.NA)

    if type(row["Precip"]) is str:
        noPrecipFormula = "Vmax ~ C(Vegetation) + C(Vegetation):C(Precip)"
        noPrecipModel = mixedlm(noPrecipFormula, df, randomFormula2,
                                groups="ID")
        noPrecipResults = noPrecipModel.fit(method=["lbfgs"])
        noPrecipResiduals = noPrecipResults.resid
        enzymeResults.loc[index, "noPrecipSSR"] = np.sum(noPrecipResiduals**2)
        modelFits.append(noPrecipResults)
    else:
        modelFits.append(pd.NA)

    if type(row["interaction"]) is str:
        noInteractionFormula = "Vmax ~ C(Vegetation) + C(Precip)"
        noInteractionModel = mixedlm(noInteractionFormula, df, randomFormula2,
                                     groups="ID")
        noInteractionResults = noInteractionModel.fit(method=["lbfgs"])
        noInteractionResiduals = noInteractionResults.resid
        enzymeResults.loc[index, "noInteractionSSR"] = np.sum(noInteractionResiduals**2)
        modelFits.append(noInteractionResults)
    else:
        modelFits.append(pd.NA)
    modelsNoSignificantFixed.append(modelFits)

# Calculating R2 for models with significant fixed effects
enzymeResults["noVegR2"] = 1 - (enzymeResults.noVegSSR/enzymeResults.nullModelSSR)
enzymeResults["noPrecipR2"] = 1 - (enzymeResults.noPrecipSSR/enzymeResults.nullModelSSR)
enzymeResults["noInteractionR2"] = 1 - (enzymeResults.noInteractionSSR/enzymeResults.nullModelSSR)

# Calculating f2 effect sizes for significant fixed effects
enzymeResults["vegf2"] = (enzymeResults.fullModelR2 - enzymeResults.noVegR2)/(1 - enzymeResults.fullModelR2)
enzymeResults["precipf2"] = (enzymeResults.fullModelR2 - enzymeResults.noPrecipR2)/(1 - enzymeResults.fullModelR2)
enzymeResults["interactionf2"] = (enzymeResults.fullModelR2 - enzymeResults.noInteractionR2)/(1 - enzymeResults.fullModelR2)
"""These effect sizes are practically zero except for the interaction effect
on PPO. What the fuck? You know what, these effect size calculations are
worthless. Let's not put them into the paper or the paper's appendix.

Instead, let's calculate Cohen's D for vegetation and precipitation if they're
significant. I can't do anything with the significant interaction for PPO.
"""


def CohenD(series1, series2):
    """
    Calculates Cohen's D for 2 different groups. Intended to be used for
    calculating effect sizes for significant vegetation and precipitation
    effects but not significant interactions of these 2 independent variables.

    Parameters
    ----------
    series1 : Pandas series
        Dependent variable values of the first group.
    series2 : Pandas series
        Dependent variable values of the second group.

    Returns
    -------
    Cohen's D representing an effect size of the difference between 2 groups.

    """
    mean1 = series1.mean()
    var1 = series1.var()
    size1 = series1.shape[0]

    mean2 = series2.mean()
    var2 = series2.var()
    size2 = series2.shape[0]

    pooledVar = (((size1 - 1)*var1) + ((size2 - 1)*var2))/(size1 + size2 - 2)
    pooledSD = pooledVar**(1/2)
    d = np.abs((mean1 - mean2)/pooledSD)
    return d


dependentVarDict = {"Enzyme": "Vmax", "functionalGroup": "spectralArea",
                    "substrate": "genesPerMillionReads",
                    "compositionParameter": "value"}


def mainEffectsCohenD(significanceLists, data, dataset):
    """
    Abstracts the calculation of Cohen's D for the main effects of each
    dataset, so that

    Parameters
    ----------
    significanceLists : List of lists
        Each list within the overall list contains the linear mixed effect
        model testing results of the independent variables and their
        interaction. All the lists describe dependent variables within a
        specific dataset (e.g. Vmax for all enzymes, spectral area for specific
        functional groups, etc).
    data : Pandas dataframe
        A specific dataset.
    dataset : string
        Conceptually, this string describes the dataset that is the 'data'
        dataframe. More specifically, this describes the specific column within
        the dataset that will be filtered for. For now, the valid values are
        "Enzyme" and "functionalGroup".

    Returns
    -------
    A Pandas dataframe containing Cohen's D for the main effects.

    """
    dependentVarCol = dependentVarDict[dataset]

    for i in range(len(significanceLists)):
        sublist = significanceLists[i]
        dependent = sublist[0]
        transformation = sublist[-1]
        queryString = "{} == @dependent".format(dataset)
        df = data.query(queryString)
        if transformation == "log10":
            df[dependentVarCol] = np.log10(df[dependentVarCol])
        elif transformation == "reciprocal":
            df[dependentVarCol] = 1/df[dependentVarCol]

        if isinstance(sublist[1], str):
            CSS = df.loc[df.Vegetation == "CSS", dependentVarCol]
            grassland = df.loc[df.Vegetation == "Grassland", dependentVarCol]

            significanceLists[i][1] = CohenD(CSS, grassland)
        if isinstance(sublist[2], str):
            reduced = df.loc[df.Precip == "Reduced", dependentVarCol]
            ambient = df.loc[df.Precip == "Ambient", dependentVarCol]

            significanceLists[i][2] = CohenD(reduced, ambient)
    results = (pd.DataFrame(significanceLists, columns=resultsDFcols)
               .astype({"interaction": "string", "Vegetation": "float64",
                        "Precip": "float64"})
               )
    results.loc[results.Vegetation.isna(), "Vegetation"] = pd.NA
    results.loc[results.Precip.isna(), "Precip"] = pd.NA
    results.loc[results.interaction.isna(), "interaction"] = pd.NA
    results = results[["dependent", "transformation", "Vegetation", "Precip",
                       "interaction"]]
    return results


enzymeResults = mainEffectsCohenD(enzymeResultsLists, VmaxDF, "Enzyme")
"""For the other datasets, let's also calculate Cohen's D for precipitation
and vegetation instead of calculating f2, and let's not calculate Cohen's D
for significant interactions. This is because Cohen's D can only be calculated
for binary categorical variables, and the interactions are 4 categories."""
# %%
# Wrangling the litter chemistry data prior to analysis

splitted = litterChemDF.id.str.split("-", expand=True)
splitted.columns = ["originalTimepoint", "originalPlot"]
splitted["originalTimepoint"] = splitted["originalTimepoint"].str.lstrip("Tt")
splitted = (splitted.astype({"originalTimepoint": "int64",
                             "originalPlot": "category"})
            )
litterChemDF["originalTimepoint"] = splitted["originalTimepoint"]
litterChemDF["ID"] = splitted["originalPlot"]

# Doing a sanity check to make sure that the time points between the original
# plot IDs and the new time points match
unmatchedTimepoints = litterChemDF.query("originalTimepoint != timePoint")
"Empty dataframe, so the original and new time points match"

litterChemDF = (litterChemDF.drop(columns="originalTimepoint")
                .rename(columns={"id": "plotTimepoint"})
                )

# Now the litter chemistry data is ready for analysis
# %%
# Creating linear mixed effect models for glycosidic bond spectral area

# Creating dataframes of spectral area for each band assignment
glycosidicBondDF = litterChemDF.query("functionalGroup == 'glycosidicBond'")
C_O_stretchDF = litterChemDF.query("functionalGroup == 'C_O_stretching'")
alkaneDF = litterChemDF.query("functionalGroup == 'alkane'")
lipidDF = litterChemDF.query("functionalGroup == 'lipid'")
carboEster1df = litterChemDF.query("functionalGroup == 'carboEster1'")
carboEster2df = litterChemDF.query("functionalGroup == 'carboEster2'")
amide1df = litterChemDF.query("functionalGroup == 'amide1'")
amide2df = litterChemDF.query("functionalGroup == 'amide2'")

# Starting the process of running models for glycosidic bond spectral area
glycosidicBond1 = mixedLinearModel(glycosidicBondDF, "spectralArea")
"""Shapiro-Wilk p = 0.50877, it's normal, so let's just use the untransformed
data."""


"""
Significant vegetation and interaction, no main precipitation effects.
"""
glycosidicBondSignificance = ["glycosidicBond", "***", None, "*", None]
# %%
# Creating models for C-O carbonyl stretches
C_O_stretch1 = mixedLinearModel(C_O_stretchDF, "spectralArea")
"Shapiro-Wilk p = 0.425"


"""
Significant vegetation effect, no precipitation or interaction effect
"""
C_O_stretchSignificance = ["C_O_stretching", "***", None, None, None]
# %%
# Creating models for alkane band assignments, which really should be lignin
alkane1 = mixedLinearModel(alkaneDF, "spectralArea")
"Shapiro-Wilk p = 0.00012165"

# Log 10 transforming and rerunning
alkane2 = mixedLinearModel(alkaneDF, "spectralArea", "log10")
"Shapiro-Wilk p = 0.0458"

# Reciprocal transforming and rerunning
alkane3 = mixedLinearModel(alkaneDF, "spectralArea", "reciprocal")
"Shapiro-Wilk p = 0.490"

"""Reciprocal transformation results in the model with the most normal
residuals. From this transformation, we get a significant vegetation effect,
and no main precipitation effect or interaction"""

alkaneSignificance = ["alkane", "*", None, None, "reciprocal"]
# %%
# Creating models for the lipid band assignment
lipid1 = mixedLinearModel(lipidDF, "spectralArea")
"Shapiro-Wilk p = 0.307"

"""
Significant vegetation effect, no interaction or precipitation effect
"""
lipidSignificance = ["lipid", "***", None, None, None]
# %%
# Creating models for carboEster1 band assignment
carboEster1model1 = mixedLinearModel(carboEster1df, "spectralArea")
"Shapiro-Wilk p = 0.0141"

# Log 10 transforming and rerunning model
carboEster1model2 = mixedLinearModel(carboEster1df, "spectralArea", "log10")
"Shapiro-Wilk p = 0.0170"

# Reciprocal transforming and rerunning
carboEster1model3 = mixedLinearModel(carboEster1df, "spectralArea",
                                     "reciprocal")
"Shapiro-Wilk p = 0.0942, let's use this"

"""From the reciprocal transformation, there is a significant vegetation
effect but no main precipitation effect or interaction"""
carboEster1significance = ["carboEster1", "***", None, None, "reciprocal"]
# %%
# Creating models for carboEster2 band assignment
carboEster2model1 = mixedLinearModel(carboEster2df, "spectralArea")
"Shapiro-Wilk p = 0.591, no need for transformation"

"""Significant vegetation effect, significant interaction"""
carboEster2significance = ["carboEster2", "***", None, "*", None]
# %%
# Creating models for the amide1 band assignment
amide1model1 = mixedLinearModel(amide1df, "spectralArea")
"Shapiro-Wilk p = 0.0963"

"""No significance at all"""
amide1significance = ["amide1", None, None, None, None]
# %%
# Creating models for the amide2 band assignment
amide2model1 = mixedLinearModel(amide2df, "spectralArea")
"Shapiro-Wilk p = 0.0908"

"""Significant vegetation effect, slightly significant precipitation effect,
slightly significant interaction."""
amide2significance = ["amide2", "*", "-", "-", None]
# %%
# Calculating effect sizes for litter chemistry functional groups
litterChemSignificance = [glycosidicBondSignificance, C_O_stretchSignificance,
                          alkaneSignificance, lipidSignificance,
                          carboEster1significance, carboEster2significance,
                          amide1significance, amide2significance]
litterChemResults = mainEffectsCohenD(litterChemSignificance, litterChemDF,
                                      "functionalGroup")
# %%
"""Ashish had given me a new metagenomic dataset in which the dependent
variable is the number of genes per million reads that degrade a specific
substrate. Not entirely sure what the units mean, still waiting for some
clarification from him. But this is what we're using now instead of the
weird previous dataset of the proportion of CAZyme genes that degrade a
specific substrate."""
metagenomicFolder = repository/"CAZyme metagenomic data"
metagenomicFileName = "CAZyme gene counts.csv"
metagenomicFilePath = metagenomicFolder/metagenomicFileName
metagenomic = (pd.read_csv(metagenomicFilePath, dtype={"timePoint": "int8"})
               .merge(timeDF, "inner", "timePoint")
               )
putativeSubstrates = metagenomic.substrate.drop_duplicates().tolist()
# %%
# Running linear mixed effect models for cellulose
celluloseDF = metagenomic.query("substrate == 'cellulose'")

# Running the model without transformations
cellulose1 = mixedLinearModel(celluloseDF, "genesPerMillionReads")

"""Both main effects are significant. Interaction is significant. Residuals are
normally distributed."""
celluloseSignificance = ["cellulose", "***", "**", "***", None]
# %%
# Running linear mixed effect models for chitin

chitinDF = metagenomic.query("substrate == 'chitin'")

# Running the model without transformations
chitin1 = mixedLinearModel(chitinDF, "genesPerMillionReads")

"""Residuals are normally distributed. Significant Vegetation effect only"""
chitinSignificance = ["chitin", "*", None, None, None]
# %%
# Running linear mixed effect models for hemicellulose
hemicelluloseDF = metagenomic.query("substrate == 'hemicellulose'")

# Running model without transformations
hemicellulose1 = mixedLinearModel(hemicelluloseDF, "genesPerMillionReads")

"Residuals are normally distributed. Significant main effects. No interaction"
hemicelluloseSignificance = ["hemicellulose", "***", "*", None, None]
# %%
# Running linear mixed effect models for lignin
ligninDF = metagenomic.query("substrate == 'lignin'")

# Running model without transformations
lignin1 = mixedLinearModel(ligninDF, "genesPerMillionReads")

"Normally distributed residuals. Significant main effects. No interaction"
ligninSignificance = ["lignin", "***", "**", None, None]
# %%
# Running linear mixed effect models for oligosaccharides
oligosaccharidesDF = metagenomic.query("substrate == 'oligosaccharides'")

# Running model without transformations
oligosaccharides1 = mixedLinearModel(oligosaccharidesDF,
                                     "genesPerMillionReads")

# Residuals not normally distributed. Transforming and re-testing
oligosaccharides2 = mixedLinearModel(oligosaccharidesDF,
                                     "genesPerMillionReads", "log10")

# Log-10 transformation still not normal. Reciprocal transforming and re-test
oligosaccharides3 = mixedLinearModel(oligosaccharidesDF,
                                     "genesPerMillionReads", "reciprocal")

"""Reciprocal transformation results in normal residuals. Marginally
significant vegetation effect. No other significance"""
oligosaccharidesSignificance = ["oligosaccharides", "-", None, None,
                                "reciprocal"]
# %%
# Running linear mixed effect models for peptidoglycan
peptidoglycanDF = metagenomic.query("substrate == 'peptidoglycan'")

# Running model without transformations
peptidoglycan1 = mixedLinearModel(peptidoglycanDF, "genesPerMillionReads")

"Normally distributed residuals. No significant main effects or interaction"
peptidoglycanSignificance = ["peptidoglycan", None, None, None, None]
# %%
# Running linear mixed effect models for polysaccharides
polysaccharidesDF = metagenomic.query("substrate == 'polysaccharides'")

# Running model without transformations
polysaccharides1 = mixedLinearModel(polysaccharidesDF, "genesPerMillionReads")

"Normally distributed residuals. No significant main effects or interaction"
polysaccharidesSignificance = ["polysaccharides", None, None, None, None]
# %%
# Running linear mixed effect models for starch
starchDF = metagenomic.query("substrate == 'starch'")

# Running model without transformations
starch1 = mixedLinearModel(starchDF, "genesPerMillionReads")

"""Residuals are, almost not normal (Shapiro-Wilk p = 0.0600993). Hmm. Let's
keep it. Significant vegetation effect. No precipitation effect. No
interaction."""
starchSignificance = ["starch", "*", None, None, None]
# %%
"""Calculating Cohen's D as effect size for significant or marginally
significant main effects"""
metagenomicSignificance = [celluloseSignificance, chitinSignificance,
                           hemicelluloseSignificance, ligninSignificance,
                           oligosaccharidesSignificance,
                           peptidoglycanSignificance,
                           polysaccharidesSignificance, starchSignificance]
metagenomicResults = mainEffectsCohenD(metagenomicSignificance, metagenomic,
                                       "substrate")
# %%
# Exporting the linear mixed model effect results
mixedModelEffectsName = "Linear mixed models, Cohen's D for main effects.xlsx"
note = ["The Vmax and FTIR data were previously analyzed using 3-way ANOVAs",
        "followed by a Tukey post-hoc. This time, they were re-analyzed using",
        "linear mixed effect models with time and plot as random effects and",
        "precipitation, vegetation, and their interaction as fixed effects.",
        "In addition, the residuals of the linear mixed effect models were",
        "tested for normality by applying the SHapiro-Wilk test on a model's",
        "residuals, while for the ANOVA analysis, the residuals were not",
        "tested for normality, and all the Vmax were log10 transformed.",
        "Another difference is that effect sizes for main effects were",
        "calculated only after performing a linear mixed effect model and",
        "before running any Tukey post-hoc while previously, they were only",
        "calculated after running Tukey post-hoc. A final difference is that",
        "Effect sizes this time were calculated at a significance level of",
        "alpha = 0.10. Tukey's post-hoc should be stricter, I'll drop the",
        "significance level down to 0.05 later.",
        "",
        "Ashish had also re-worked his metagenomic dataset to sum up the",
        "number of CAZyme genes dedicated to a specific substrate per million",
        "reads, as opposed to last time in which he calculated a proportion",
        "of CAZyme genes dedicated to a specific substrate out of all",
        "possible CAZyme genes. I analyzed this dataset using linear mixed",
        "effect models as well."]
readMe = pd.DataFrame({"Note": note})
statsFolder = repository/"Statistical analyses"
resultsPath = statsFolder/mixedModelEffectsName
if os.path.exists(resultsPath) is False:
    with pd.ExcelWriter(resultsPath) as w:
        readMe.to_excel(w, "ReadMe", index=False)
        enzymeResults.to_excel(w, "Vmax", index=False)
        litterChemResults.to_excel(w, "litterChem", index=False)
        metagenomicResults.to_excel(w, "CAZyme metagenomic", index=False)
# %%
"""Analyzing community composition data from Ashish."""
compositionFilePath = metagenomicFolder/"Community composition.xlsx"
compositionData = (pd.read_excel(compositionFilePath, "Data")
                   .merge(timeDF, "inner", "timePoint")
                   )
bacterialReads = compositionData.query("compositionParameter == 'bac'")
fungalReads = compositionData.query("compositionParameter == 'fun'")
fb = compositionData.query("compositionParameter == 'fb'")
shannon = compositionData.query("compositionParameter == 'shan'")
# %%
bacterialReads1 = mixedLinearModel(bacterialReads, "value")
# Shapiro-Wilk p = 0.001329

bacterialReads2 = mixedLinearModel(bacterialReads, "value", "log10")
# Shapiro-Wilk p = 6.4352*10^-5

bacterialReads3 = mixedLinearModel(bacterialReads, "value", "square root")
# Shapiro-Wilk p = 0.00028986

bacterialReads3 = mixedLinearModel(bacterialReads, "value", "reciprocal")
# Shapiro-Wilk p = 2.793*10^-6

"""So the untransformed is the most normal. Vegetation is marginally
significant while precipitation and the interaction aren't."""
bacterialSignificance = ["bac", "-", None, None, None]
# %%
fungalReads1 = mixedLinearModel(fungalReads, "value")
# Shapiro-Wilk p = 0.001009

fungalReads2 = mixedLinearModel(fungalReads, "value", "log10")
# Shapiro-Wilk p = 0.3388

fungalReads3 = mixedLinearModel(fungalReads, "value", "square root")
# Shapiro-Wilk p = 0.2484

fungalReads4 = mixedLinearModel(fungalReads, "value", "reciprocal")
# Shapiro-Wilk p = 0.00206

"""Let's go with the log-10 transformation. p < 0.05 for vegetation, no
significance for precipitation and the interaction"""
fungalSignificance = ["fun", "*", None, None, "log10"]
# %%
fb1 = mixedLinearModel(fb, "value")
# Shapiro-Wilk p = 2.676*10^-6

fb2 = mixedLinearModel(fb, "value", "log10")
# Shapiro-Wilk p = 0.5429

fb3 = mixedLinearModel(fb, "value", "square root")
# Shapiro-Wilk p = 0.000719

fb4 = mixedLinearModel(fb, "value", "reciprocal")
# Shapiro-Wilk p = 0.001969

"""With the log-10 transformation, it's, quite weird that there's no
p-value estimate for vegetation. No significance for precipitation and the
interaction. Vegetation is marginally significant under a square root
transformation and is < 0.001 under a reciprocal transformation. Guess I'll go
with the reciprocal transformation for now."""
fbSignificance = ["fb", "***", None, None, "reciprocal"]
# %%
shannon1 = mixedLinearModel(shannon, "value")
# Shapiro-Wilk p = 2.411*10^-6

shannon2 = mixedLinearModel(shannon, "value", "log10")
# Shapiro-Wilk p = 1.040*10^-6

shannon3 = mixedLinearModel(shannon, "value", "square root")
# Shapiro-Wilk p = 1.548*10^-6

shannon4 = mixedLinearModel(shannon, "value", "reciprocal")
# shapiro-Wilk p = 3.2134*10^-7

"""Going with the untransformed data. Vegetation < 0.001, no significance for
precipitation, p = 0.104 for interaction. I'll do Tukey's on vegetation and
the interaction."""
shannonSignificance = ["shan", "***", None, "-", None]
# %%
# Calculating Cohen's D for community composition
compositionSignificance = [bacterialSignificance, fungalSignificance,
                           fbSignificance, shannonSignificance]
compositionResults = mainEffectsCohenD(compositionSignificance, compositionData, "compositionParameter")
# %%
"""Exporting the community composition linear mixed modeling results and
Cohen's D"""
compositionResultsPath = statsFolder/"Linear mixed models, community composition.xlsx"
note2 = ["This file contains results from mixed effect modeling of 'bac'",
         "(which I believe refers to bacterial reads), 'fun' (which likely",
         "refers to fungal reads), 'fb' (fungal:bacterial ratios calculated",
         "by dividing 'fun' by 'bac'), and 'shan' (Shannon's taxonomic",
         "diversity). I calculated Cohen's D for significant or marginally",
         "significant main effects, marginally significant being when p is",
         "between 0.05 and 0.10 or is just slightly above 0.10. I would have",
         "used asterisks to describe significant interactions, although there",
         "weren't any significant interactions"]
readMe2 = pd.DataFrame({"note": note2})
if os.path.exists(compositionResultsPath) is False:
    with pd.ExcelWriter(compositionResultsPath) as w:
        readMe2.to_excel(w, "ReadMe", index=False)
        compositionResults.to_excel(w, "Community composition", index=False)
