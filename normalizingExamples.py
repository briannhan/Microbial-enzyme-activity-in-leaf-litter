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
py.style.use("dark_background")
# %%
# Equilibrium Chemistry Approximation equation


def ECA(k2, E, Km, substrateConcentration):
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
totalActivity = ECA(k2, enzyme, Km, substrate)  # mol L^-1 s^-1
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
activityModk2 = ECA(normk2, enzyme, Km, substrate)  # mol g^-1 L^-1 s^-1
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
activityModSEKm = ECA(k2, normEnzyme, normKm, normSubstrate)
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
normActivity3Modk2 = ECA(normk2Sample3, enzyme, Km, substrate)
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
normActivity3ModSEKm = ECA(k2, normEnzyme3, normKm3, normSubstrate3)
py.subplot(3, 3, 9)
py.xlabel("Normalized substrate concentration (mole g^-1 L^-1)")
py.title("3, Normalized enzyme activity & substrates (S, E, & Km/litterMass)")
py.ylim(top=25)
py.xlim(right=750)
py.plot(normSubstrate3, normActivity3ModSEKm)

figure2Path = figuresFolder/(figure2Title + ".png")
py.savefig(figure2Path)
