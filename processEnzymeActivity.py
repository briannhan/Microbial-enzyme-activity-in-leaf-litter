# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:18:28 2021

@author: Brian Chung
The purpose of this script is to process enzyme activity that had been
calculated using the calculateEnzymeActivity.py script. This script uses an
Excel spreadsheet file I've made containing my labels for how I processed the
errors I saw when looking at the unprocessed enzyme activity graphs.

I envision that this script will process enzyme activity in a series of steps.
After each step, the data will be fitted to the single-substrate,
single-enzyme formulation of an approximation of equilibrium chemistry (ECA)
according to Tang and Riley 2013. The fitting of the data will be constrained
so that values of k2 (the maximum product genesis rate), Km (Michaelis-Menten
constant), and enzyme concentration are all positive.
"""
