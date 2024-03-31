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
metagenomic.Precip = metagenomic.Precip.str.strip()
metagenomic.Vegetation = metagenomic.Vegetation.str.strip()
metagenomic.id = metagenomic.id.str.strip()

metadataFolder = repository/"Metadata"
metadataFileName = "2020-11-10-Litter-bag-codes.xlsx"
metadataPath = metadataFolder/metadataFileName
metadataNames = ["Vegetation", "Precip", "timePoint", "ID"]
metadata = (pd.read_excel(metadataPath, "Wet up expt", names=metadataNames,
                          dtype="string")
            .drop(columns="timePoint")
            .drop_duplicates()
            )
metadata.loc[metadata.Vegetation.str.contains("Grass"), "Vegetation"] = "Grassland"
metadata.loc[metadata.Vegetation.str.contains("Shrub"), "Vegetation"] = "CSS"
metadata.Precip = metadata.Precip.str.strip()
# %%
metagenomic["idDerivedTimePoint"] = metagenomic.id.str[1]

metagenomic["ID"] = metagenomic.id.str[2:]


metagenomic["sampDerivedTimePoint"] = metagenomic.samp.str[1]
metagenomic.loc[metagenomic.sampDerivedTimePoint == "1", "sampDerivedTimePoint"] = "0"
metagenomic.loc[metagenomic.sampDerivedTimePoint == "3", "sampDerivedTimePoint"] = "5"
metagenomic.loc[metagenomic.sampDerivedTimePoint == "2", "sampDerivedTimePoint"] = "3"
metagenomic.loc[metagenomic.sampDerivedTimePoint == "4", "sampDerivedTimePoint"] = "6"

# Seeing if the 2 derived timepoints don't agree with each other
unequalTimepoints = metagenomic.query("idDerivedTimePoint != sampDerivedTimePoint")
"""Purposefully switched the order of the lines that assigned 3 to 5 and 2 to
3 around because if I didn't switch the order, then there would be an extra
16 plots assigned to timepoint 5 and timepoint 3 wouldn't have any plots
left.

Anyway, the time points are all in agreement.
"""
metagenomic = (metagenomic.drop(columns=["id", "sampDerivedTimePoint"])
               .rename(columns={"idDerivedTimePoint": "timePoint"})
               )

# Checking if the vegetation and precipitation columns are accurate
metagenomic.loc[metagenomic.Vegetation == "grass", "Vegetation"] = "Grassland"
metagenomic.loc[metagenomic.Vegetation == "shrub", "Vegetation"] = "CSS"
metagenomic.loc[metagenomic.Precip == "amb", "Precip"] = "Ambient"
metagenomic.loc[metagenomic.Precip == "red", "Precip"] = "Reduced"
independentMergeTest = (metagenomic.merge(metadata, "outer",
                                          ["ID", "Vegetation", "Precip"],
                                          indicator=True)
                        .query("_merge != 'both'")
                        )
"""Vegetation and precipitation columns are also accurate"""

# %%
metagenomic = metagenomic.melt(["ID", "timePoint", "samp", "Vegetation",
                                "Precip"], var_name="substrate",
                               value_name="genesPerMillionReads")
exportName = "CAZyme gene counts.csv"
exportPath = metagenomicFolder/exportName
if os.path.exists(exportName) is False:
    metagenomic.to_csv(exportPath, index=False)
