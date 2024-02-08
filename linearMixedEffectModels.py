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
litterChemDtype = {"id": "string", "Vegetation": "string", "Precip": "string",
                   "timePoint": "category", "functionalGroup": "string"}
litterChemDF = (pd.read_excel(litterChemPath, 0, dtype=litterChemDtype)
                .drop_duplicates()
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
# timeDict = {"timePoint": ["deployment", "0", "3", "5", "6"]}
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
# AG.timePoint = AG.timePoint.astype("int64")
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
resultsDFcols = ["enzyme", "vegetation", "precipitation", "twoWay"]
AGsignificance = ["AG", None, None, None]
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

APsignificance = ["AP", None, "-", None]
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

    model = mixedlm(formula, data, randomFormula2, groups="ID")
    modelResults = model.fit(method=["lbfgs"])
    modelSummary = modelResults.summary()
    modelResiduals = modelResults.resid
    modelNormality = stats.shapiro(modelResiduals)
    return modelSummary, modelNormality


# %%
# Creating mixed linear models for beta-glucosidase, BG
BG = VmaxDF.query("Enzyme == 'BG'")
BGmodel1summary, BGmodel1normality = mixedLinearModel(BG, "Vmax")
"""Shapiro-Wilk p-value after converting time = 4.527*10^-16"""

# Log 10 transforming BG and re-testing
BGmodel2summary, BGmodel2normality = mixedLinearModel(BG, "Vmax", "log10")
"""Shapiro-Wilk p-value after converting time = 1.3198*10^-10"""


# Square root transforming BG and re-testing
BGmodel3summary, BGmodel3normality = mixedLinearModel(BG, "Vmax",
                                                      "square root")
"""Square root transformation has better normality than log10, although still
not quite normal

After converting time, Shapiro-Wilk p-value = 6.1116*10^-15
"""

# Reciprocal transforming BG and re-testing
BGmodel4summary, BGmodel4normality = mixedLinearModel(BG, "Vmax", "reciprocal")
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
BGsignificance = ["BG", "***", "*", None]
# %%
# Creating mixed linear models for beta-xylosidase, BX
BX = VmaxDF.query("Enzyme == 'BX'")
BXmodel1summary, BXmodel1normality = mixedLinearModel(BX, "Vmax")
"""Convert time -> numerical, SW p-value = 3.829*10^-14"""

# Log 10 transforming BX and re-testing
BXmodel2summary, BXmodel2normality = mixedLinearModel(BX, "Vmax", "log10")
"""Improves normality, but still far from normal

Convert time -> numerical, SW p-value = 1.303*10^-05
"""

# Square root transforming BX and re-testing
BXmodel3summary, BXmodel3normality = mixedLinearModel(BX, "Vmax",
                                                      "square root")
"""Improves normality but worse than log10

Convert time -> numerical, SW p-value = 2.0909*10^-10
"""

# Reciprocal transforming BX and re-testing
BXmodel4summary, BXmodel4normality = mixedLinearModel(BX, "Vmax", "reciprocal")
"""Actually worsens normality

Convert time -> numerical, SW p-value = 6.72*10^-13
"""

"""Let's just use results from log 10 transformation.
No effects by either vegetation or precipitation, no interaction either"""
BXsignificance = ["BX", None, None, None]
# %%
# Creating mixed linear models for cellobiohydrolase, CBH
CBH = VmaxDF.query("Enzyme == 'CBH'")
CBHmodel1summary, CBHmodel1normality = mixedLinearModel(CBH, "Vmax")
"""Shapiro-Wilk p-value = 1.11*10^-6

After converting time to numerical variable, Shapiro-Wilk p-value = 0.00020449
"""

# Log 10 transforming CBH and re-testing
CBHmodel2summary, CBHmodel2normality = mixedLinearModel(CBH, "Vmax", "log10")
"""Shapiro-Wilk p-value = 0.6667, let's keep this transformation.
Significant vegetation effect, significant precipitation effect,
no interaction

Shapiro-Wilk p-value after converting time to continuous variable = 0.495
"""

"""From log10 transformation, significant vegetation effect, marginally
significant precipitation effect that warrants at least a Tukey. No significant
interaction."""
CBHsignificance = ["CBH", "***", "-", None]
# %%
# Creating mixed linear models for LAP
LAP = VmaxDF.query("Enzyme == 'LAP'")
LAPmodel1summary, LAPmodel1normality = mixedLinearModel(LAP, "Vmax")
"""Shapiro-Wilk p = 4.503*10^-12"""

# Log 10 transforming and re-testing
LAPmodel2summary, LAPmodel2normality = mixedLinearModel(LAP, "Vmax", "log10")
"""Shapiro-Wilk p = 6.254*10^-5"""

# Square root transforming and re-testing
LAPmodel3summary, LAPmodel3normality = mixedLinearModel(LAP, "Vmax",
                                                        "square root")
"""Shapiro-Wilk p = 1.81*10^-9"""

# Reciprocal transforming and re-testing
LAPmodel4summary, LAPmodel4normality = mixedLinearModel(LAP, "Vmax",
                                                        "reciprocal")
"""Shapiro-Wilk p = 0.000354, the most normal transformation. Let's use
this"""

"""From reciprocal transformation, significant vegetation effect, no
precipitation effect, no interaction"""
LAPsignificance = ["LAP", "**", None, None]
# %%
# Creating mixed linear models for NAG
NAG = VmaxDF.query("Enzyme == 'NAG'")
NAGmodel1summary, NAGmodel1normality = mixedLinearModel(NAG, "Vmax")
"""Shapiro-Wilk normality p = 5.189*10^-7"""

# Log 10 transforming and re-testing
NAGmodel2summary, NAGmodel2normality = mixedLinearModel(NAG, "Vmax", "log10")
"""Shapiro-Wilk p = 0.000451"""

# Square root transforming and re-testing
NAGmodel3summary, NAGmodel3normality = mixedLinearModel(NAG, "Vmax",
                                                        "square root")
"""Shapiro-Wilk normality p = 2.778*10^-6"""

# Reciprocal transforming and re-testing
NAGmodel4summary, NAGmodel4normality = mixedLinearModel(NAG, "Vmax",
                                                        "reciprocal")
"""Normality worsens, SW p = 1.0003*10^-9"""

"""From log 10 transformation, significant Vegetation effect, no precipitation
effect, no interaction"""
NAGsignificance = ["NAG", "***", None, None]
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
PPOsignificance = ["PPO", None, None, "**"]

"""So now all the enzymes have been tested. Moving onto the litter chemistry
functional groups."""
# %%