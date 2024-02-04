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
                   "timePoint": "string", "functionalGroup": "string"}
litterChemDF = (pd.read_excel(litterChemPath, 0, dtype=litterChemDtype)
                .drop_duplicates()
                )
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
randomFormula = "~0 + C(timePoint)"
randomFormula2 = "~C(timePoint)"
randomFormula3 = "~C(timePoint) - 1"

"""Formulas 3 and 1 remove the intercept for timePoint and only estimated the
slope for timePoint. Which, might not be that good of an idea"""
AGmodel1 = mixedlm(fullFormula, AG, randomFormula, groups="ID")
AGmodel1results = AGmodel1.fit(method=["lbfgs"])
AGmodel1summary = AGmodel1results.summary()

AGmodel2 = mixedlm(fullFormula, AG, randomFormula2, groups="ID")
AGmodel2results = AGmodel2.fit(method=["lbfgs"])
AGmodel2summary = AGmodel2results.summary()

AGmodel3 = mixedlm(fullFormula, AG, randomFormula3, groups="ID")
AGmodel3results = AGmodel3.fit(method=["lbfgs"])
AGmodel3summary = AGmodel3results.summary()

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
"""Model residuals still not normal after log10 transformation"""

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
precipitation"""

"""Checking model residuals to see if they're likely from a normal population
"""
APmodel1residuals = APmodel1results.resid
APmodel1normality = stats.shapiro(APmodel1residuals)
"""P-value of Shapiro-Wilk is 1.5*10^-14, very small, indicating it's unlikely
that the residuals are from a normal population. Time to transform it."""

# Log10 transforming AP and running model again
APlog10 = AP.copy()
APlog10["Vmax"] = np.log10(APlog10["Vmax"])
APmodel2 = mixedlm(fullFormula, APlog10, randomFormula2, groups="ID")
APmodel2results = APmodel2.fit(method=["lbfgs"])
APmodel2summary = APmodel2results.summary()
APmodel2residuals = APmodel2results.resid
APmodel2normality = stats.shapiro(APmodel2residuals)
"""Still not quite normal, with Shapiro-Wilk normality p-value = 0.000394"""

# Natural log transform and running model again
APln = AP.copy()
APln["Vmax"] = np.log10(APln["Vmax"])
APmodel3 = mixedlm(fullFormula, APln, randomFormula2, groups="ID")
APmodel3results = APmodel3.fit(method=["lbfgs"])
APmodel3summary = APmodel3results.summary()
APmodel3residuals = APmodel3results.resid
APmodel3normality = stats.shapiro(APmodel3residuals)
"""Normality identical to log10 transformations. I'm not gonna do natural log
transformations again in the future, then."""

# Square root transform and running again
APsquareRoot = AP.copy()
APsquareRoot["Vmax"] = APsquareRoot["Vmax"]**(1/2)
APmodel4 = mixedlm(fullFormula, APsquareRoot, randomFormula2, groups="ID")
APmodel4results = APmodel4.fit(method=["lbfgs"])
APmodel4summary = APmodel4results.summary()
APmodel4residuals = APmodel4results.resid
APmodel4normality = stats.shapiro(APmodel4residuals)
"""Normality improves slightly but still worse than log transformations"""

# Reciprocal transform and running again
APreciprocal = AP.copy()
APreciprocal["Vmax"] = 1/APreciprocal["Vmax"]
APmodel5 = mixedlm(fullFormula, APreciprocal, randomFormula2, groups="ID")
APmodel5results = APmodel5.fit(method=["lbfgs"])
APmodel5summary = APmodel5results.summary()
APmodel5residuals = APmodel5results.resid
APmodel5normality = stats.shapiro(APmodel5residuals)
"""Normality improves, p-value = 0.0002177, lower p-value than log
transformations, so still slightly lower normality than log transformations"""

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
disappeared."""

APsignificance = ["AP", None, None, None]
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

# Log 10 transforming BG and re-testing
BGmodel2summary, BGmodel2normality = mixedLinearModel(BG, "Vmax", "log10")

# Square root transforming BG and re-testing
BGmodel3summary, BGmodel3normality = mixedLinearModel(BG, "Vmax",
                                                      "square root")
"""Square root transformation has better normality than log10, although still
not quite normal"""

# Reciprocal transforming BG and re-testing
BGmodel4summary, BGmodel4normality = mixedLinearModel(BG, "Vmax", "reciprocal")
"""Improves normality but worse than log10 and square root. Let's use results
from square root transformation"""

"""From square root transformation, effects by vegetation and precipitation
disappeared."""
BGsignificance = ["BG", None, None, None]
# %%
# Creating mixed linear models for beta-xylosidase, BX
BX = VmaxDF.query("Enzyme == 'BX'")
BXmodel1summary, BXmodel1normality = mixedLinearModel(BX, "Vmax")

# Log 10 transforming BX and re-testing
BXmodel2summary, BXmodel2normality = mixedLinearModel(BX, "Vmax", "log10")
"""Improves normality, but still far from normal"""

# Square root transforming BX and re-testing
BXmodel3summary, BXmodel3normality = mixedLinearModel(BX, "Vmax",
                                                      "square root")
"""Improves normality but worse than log10"""

# Reciprocal transforming BX and re-testing
BXmodel4summary, BXmodel4normality = mixedLinearModel(BX, "Vmax", "reciprocal")
"""Actually worsens normality"""

"""Let's just use results from log 10 transformation.
No effects by either vegetation or precipitation, no interaction either"""
BXsignificance = ["BX", None, None, None]
# %%
