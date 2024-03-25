# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 21:01:43 2024

@author: LMAOXD

This script wrangles the new metagenomic dataset Ashish gave me, in which the
dependent variable is the number of CAZyme genes that target a specific
substrate per million reads in a plot. Still not entirely sure what the units
mean, but waiting for his response.

The following changes will be made
- id column will be separated into the plot and time point
- independent variable columns will be renamed, column names & values
- substrate columns will be melted down

I will perform the following checks on the data for internal consistency
- the vegetation and precipitation columns match with the original id column
- the original time point column matches with the original id column
"""
import os
from pathlib import Path
import pandas as pd

repository = Path(os.getcwd())
metagenomicFolder = repository/"CAZyme metagenomic data"
metagenomicFileName = "cazy_tpm_substrate_brian.csv"
metagenomicPath = metagenomicFolder/metagenomicFileName
metagenomic = (pd.read_csv(metagenomicPath)
               .rename(columns={"veg": "Vegetation", "prec": "Precip"})
               )
# %%
metagenomic["timePoint"] = metagenomic.id.str[1]

metagenomic["ID"] = metagenomic.id.str[2:]


metagenomic["ogTimePoint"] = metagenomic.samp.str[1]
metagenomic.loc[metagenomic.ogTimePoint == "1", "ogTimePoint"] = "0"
metagenomic.loc[metagenomic.ogTimePoint == "2", "ogTimePoint"] = "3"
metagenomic.loc[metagenomic.ogTimePoint == "3", "ogTimePoint"] = "5"
metagenomic.loc[metagenomic.ogTimePoint == "4", "ogTimePoint"] = "6"
