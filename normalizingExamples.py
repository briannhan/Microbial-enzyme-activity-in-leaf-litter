# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:49:35 2020

@author: Brian Chung
The purpose of this script is to further my understanding of how I should work
out the units of enzyme activity when I'm doing nonlinear regression. This
script is meant to benefit me, so apologies if future readers cannot understand
it. This script makes 2 figures. The first figure will have 5 subplots and the
second will have 9. The enzyme activities follow the equilibrium chemistry
approximation relationship of enzyme activity.
"""
import matplotlib.pyplot as py
import numpy as np
from pathlib import Path
import os
from scipy.optimize import curve_fit
py.style.use("dark_background")
# %%
# Equilibrium Chemistry Approximation equation


def ECA(substrateConcentration, k2, E, Km):
    """Calculates enzyme activity based on the Equilibrium Chemistry
    Approximation

    Parameters
    ----------
    k2 : float
        Catalytic rate constant.
    E : float
        Concentration of enzymes (mol/L)
    Km : float
        Michaelis-Menten constant (mol/L)
    substrateConcentration : numpy array of floats
        Substrate concentrations (mol/L)

    Returns
    -------
    Enzyme activity : float
        Enzyme activity (mole L^-1 s^-1)

    """
    return (k2*E*substrateConcentration)/(Km + E + substrateConcentration)


# %%
# Making figure 1
substrate = np.linspace(start=0, stop=500)  # mol/L
litterMass = 1.5  # grams
k2 = 100  # s^-1
Km = 100  # mol/L
enzyme = 0.2  # mol/L
totalActivity = ECA(substrate, k2, enzyme, Km)  # mol L^-1 s^-1
figure1Title = "Normalizing 1 enzyme sample"
py.figure(num=figure1Title, figsize=(20, 15))
# Making subplot of total enzyme activity
py.subplot(2, 3, 1)
py.xlabel("Substrate concentration (mol/L)")
py.ylabel("Total enzyme activity (mole L^-1 s^-1)")
py.title("Total enzyme activity")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate, totalActivity)

# Normalizing enzyme activity
normActivity = totalActivity/litterMass  # mol g^-1 L^-1 s^-1
py.subplot(2, 3, 2)
py.xlabel("Substrate concentration (mol/L)")
py.ylabel("Normalized enzyme activity (mole g^-1 L^-1 s^-1)")
py.title("Normalized enzyme activity")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate, normActivity)

# Dividing k2 by litterMass
normk2 = k2/litterMass  # s^-1 g^-1
activityModk2 = ECA(substrate, normk2, enzyme, Km)  # mol g^-1 L^-1 s^-1
py.subplot(2, 3, 3)
py.xlabel("Substrate concentration (mol/L)")
py.ylabel("Enzyme activity (mole g^-1 L^-1 s^-1) with mod. k2")
py.title("k2/litterMass")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate, activityModk2)

# Normalizing both enzyme activity & substrate concentration
normSubstrate = substrate/litterMass  # mol g^-1 L^-1
py.subplot(2, 3, 4)
py.xlabel("Normalized substrate concentration (mole g^-1 L^-1)")
py.ylabel("Normalized enzyme activity (mole g^-1 L^-1 s^-1)")
py.title("Normalized enzyme activity & substrates")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate, normActivity)

# Normalizing S, E, & Km
normEnzyme = enzyme/litterMass  # mol g^-1 L^-1
normKm = Km/litterMass  # mol g^-1 L^-1
activityModSEKm = ECA(normSubstrate, k2, normEnzyme, normKm)
py.subplot(2, 3, 5)
py.xlabel("Normalized substrate concentration (mole g^-1 L^-1)")
py.ylabel("Enzyme activity (mole g^-1 L^-1 s^-1) with mod. S, E, Km")
py.title("S, E, & Km/litterMass")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate, activityModSEKm)

# plotting normalized substrate activity (calculated by normalizing S, E, & Km)
# against normalized substrate activity (calculated by normalizing k2) against
# each other to verify that they are the same values
py.subplot(2, 3, 6)
py.xlabel("Normalized enzyme activity (from normalizing S, E, & Km)")
py.ylabel("Normalized enzyme activity (from normalizing k2)")
py.title("Normalized enzyme activity verification")
py.plot(activityModSEKm, activityModk2)

workingDir = Path(os.getcwd())
figuresFolder = workingDir/"Normalizing examples"
figure1Path = figuresFolder/(figure1Title + ".png")
py.savefig(figure1Path)
# %%
# Making 2nd figure where I normalize 3 samples
figure2Title = "Normalizing 3 litter samples"
totalActivity2 = totalActivity
substrate2 = substrate
litterMass2 = 1.1  # grams
py.figure(num=figure2Title, figsize=(30, 20))
# Making subplot of total enzyme activity for sample 1
py.subplot(3, 3, 1)
py.ylabel("Enzyme activity (mole L^-1 s^-1)")
py.title("1, Total activity")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate, totalActivity)

# Making subplot of normalized enzyme activity for sample 1
py.subplot(3, 3, 2)
py.ylabel("Normalized enzyme activity (mole g^-1 L^-1 s^-1)")
py.title("1, Normalized enzyme activity (k2/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate, normActivity)

# Normalizing both enzyme activity and substrates of sample 1
py.subplot(3, 3, 3)
py.ylabel("Normalized enzyme activity (mole g^-1 L^-1 s^-1)")
py.title("1, Normalized enzyme activity & substrates (S, E, & Km/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate, activityModSEKm)

# Plotting total enzyme activity for sample 2
py.subplot(3, 3, 4)
py.ylabel("Enzyme activity (mole L^-1 s^-1)")
py.title("2, Total activity")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate2, totalActivity2)

# Normalizing only enzyme activity of sample 2
normActivity2 = totalActivity2/litterMass2
py.subplot(3, 3, 5)
py.title("2, Normalized enzyme activity (k2/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate2, normActivity2)

# Sample 2, normalizing both enzyme activity & substrate concentrations
normSubstrate2 = substrate2/litterMass2
py.subplot(3, 3, 6)
py.title("2, Normalized enzyme activity & substrates (S, E, & Km/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate2, normActivity2)

# Plotting total enzyme activity for sample 3
substrate3 = substrate
totalActivity3 = totalActivity
py.subplot(3, 3, 7)
py.xlabel("Substrate concentration (mole L^-1)")
py.ylabel("Enzyme activity (mole L^-1 s^-1)")
py.title("3, Total activity")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate3, totalActivity3)

# Normalizing only enzyme activity of sample 3
litterMass3 = 0.7  # grams
normk2Sample3 = k2/litterMass3
normActivity3Modk2 = ECA(substrate, normk2Sample3, enzyme, Km)
py.subplot(3, 3, 8)
py.xlabel("Substrate concentration (mol/L)")
py.title("3, Normalized enzyme activity (k2/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(substrate3, normActivity3Modk2)

# Sample 3, normalizing both enzyme activity & substrate concentrations
normSubstrate3 = substrate3/litterMass3
normEnzyme3 = enzyme/litterMass3
normKm3 = Km/litterMass3
normActivity3ModSEKm = ECA(normSubstrate3, k2, normEnzyme3, normKm3)
py.subplot(3, 3, 9)
py.xlabel("Normalized substrate concentration (mole g^-1 L^-1)")
py.title("3, Normalized enzyme activity & substrates (S, E, & Km/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate3, normActivity3ModSEKm)

figure2Path = figuresFolder/(figure2Title + ".png")
py.savefig(figure2Path)

# %%
# In this section, I will have a single litter sample with known normalized
# values of k2, enzyme concentration, and Km. I will use these normalized
# values while changing the values of the litter mass. This is to test what
# happens to k2 and Km when I perform actual assays with varying litter masses.
# I will use the normalized parameters to calculate normalized enzyme activity,
# and then multiply the normalized enzyme activity by different litter masses
# to obtain total enzyme activity at various litter masses. I will then perform
# nonlinear regression between total enzyme activity and substrate
# concentrations to obtain non-normalized k2, enzyme concentration, and Km,
# and then divide these parameters by the litter masses to see if they are the
# same as known values of normalized parameters.
litterMass1 = 0.4  # grams
litterMass2 = 0.8  # grams
litterMass3 = 1.2  # grams
k2 = 100  # s^-1
normE_og = 20  # mol L^-1 g^-1
normKm_og = 100  # mol L^-1 g^-1
normSubstrate = np.linspace(start=0, stop=500, num=100)
normActivity = ECA(normSubstrate, k2, normE_og, normKm_og)
normParams, popCov = curve_fit(ECA, normSubstrate, normActivity,
                               bounds=([90, 15, 90], [110, 25, 110]))
py.figure(num=3, figsize=(30, 20))
py.subplot(1, 2, 1)
py.plot(normSubstrate, normActivity, 'bo')
calNormActi = ECA(normSubstrate, normParams[0], normParams[1], normParams[2])
py.plot(normSubstrate, calNormActi, 'r+')

# And imma use different arrays of normalized substrates depending on the
# litter masses
# And now the math begins

# normalized litter 1
# normSubstrate1 = substrate/litterMass1
totalActivity1 = normActivity*litterMass1
substrate1 = normSubstrate*litterMass1
params1, popCov1 = curve_fit(ECA, substrate1, totalActivity1,
                             bounds=([95, 6, 35], [105, 10, 45]))
normPara1 = params1/litterMass1
# If normalizing is ok, then I expect that the k2 in params1 will be equal
# to the k2 set above, and normalized E & Km will be equal to normE_og &
# normKm_og, respectively

# Well, doesn't seem like normalizing works out fantastically.

# normalized litter 2
# normSubstrate2 = substrate/litterMass2
# normActivity2 = ECA(k2, normE_og, normKm_og, normSubstrate2)
totalActivity2 = normActivity*litterMass2
substrate2 = normSubstrate*litterMass2
params2, popCov2 = curve_fit(ECA, substrate2, totalActivity2,
                             bounds=([95, 14, 75], [105, 18, 85]))
normPara2 = params2/litterMass2

# normSubstrate3 = substrate/litterMass3
totalActivity3 = normActivity*litterMass3
substrate3 = normSubstrate*litterMass3
params3, popCov3 = curve_fit(ECA, substrate3, totalActivity3,
                             bounds=([95, 20, 115], [105, 28, 125]))
normPara3 = params3/litterMass3
'''So, in conclusion, my prediction hold true. As a recap, this section assays
a single sample 3 times with each assay having a different litter mass. The
purpose of this section is to see if normalizing parameters is a valid way to
pre-process data prior to statistical comparisons. My hypothesis for this
section is that yes, normalizing enzyme parameters is a valid way to
pre-process data. My predictions for this hypothesis are that the
non-normalized parameters would vary from assay despite all assays being from
the same litter sample and so should have no variation in enzyme parameters,
but normalized parameters would be the same across all assays. And the results
partially follow my prediction, so my hypothesis is partially validated.
k2 remains the same if it's not normalized, contrary to my hypothesis.
However, enzyme concentrations and Km remain the same if it's normalized,
following my hypothesis. So, before making statistical comparisons, enzyme
concentrations and Km must be normalized. However, k2 doesn't need to be.'''
