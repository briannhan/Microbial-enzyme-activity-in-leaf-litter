# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 01:07:16 2025

This script wrangles the community composition data that I received from
Ashish. I received what I believed are the number of fungal and bacterial reads
which Ashish used to calculate fungal:bacterial ratios. I also received
Shannon's taxonomic diversity. I'm going to wrangle the data so that I can run
linear mixed effect models on it.
"""
import os
from pathlib import Path
import pandas as pd

repository = Path(os.getcwd())
folderPath = repository/"CAZyme metagenomic data"
readsPath = folderPath/"fb.csv"
readsData = (pd.read_csv(readsPath)
             .rename(columns={"Unnamed: 0": "plotTime"})
             )
shannonPath = folderPath/"shannon.taxa.csv"
shannonData = (pd.read_csv(shannonPath)
               .rename(columns={"Unnamed: 0": "plotTime"})
               )
# %%
combined = (readsData.merge(shannonData, "outer",
                            ["plotTime", "veg", "prec", "samp"])
            .rename(columns={"veg": "Vegetation", "prec": "Precip"})
            )

"""Out of curiosity, checking to see if dividing fungal reads by bacterial
reads does give the F:B ratios that Ashish calculated"""
combined["calculatedFB"] = combined.fun/combined.bac
combined["difference"] = combined.fb - combined.calculatedFB
"""Yeah. Divbiding fungal reads by bacterial reads does give the ratios that he
calculated."""

combined = combined.drop(columns=["calculatedFB", "difference"])
combined["timePoint"] = combined.plotTime.str[1]
combined["IDandRemainder"] = combined.plotTime.str[2:]
combined = combined.astype({"timePoint": "int8"})
plotDF = (combined.IDandRemainder.str.split("_", expand=True)
          .rename(columns={0: "ID", 1: "remainder"})
          )
combined = combined.join(plotDF)

"""Not quite sure what the digit following the plot number is. I'll just delete
it. Doesn't seem important"""
combined = (combined.drop(columns=["IDandRemainder", "remainder"])
            .melt(["plotTime", "Vegetation", "Precip", "samp", "timePoint", "ID"],
                  var_name="compositionParameter")
            )

combined.loc[combined.Vegetation == "Grass", "Vegetation"] = "Grassland"
combined.loc[combined.Vegetation == "Shrub", "Vegetation"] = "CSS"
# %%
"""Exporting the wrangled data"""
note = ["This is the community composition data that Ashish used to make the",
        "figures on Shannon's taxonomic diversity and fungal:bacterial",
        "ratios. He sent the data over to me, and I wrangled the data to be",
        "in a format that allows for statistical analysis using my existing",
        "code.",
        "",
        "The data is in the spreadsheet called 'data'. The column 'plotTime'",
        "contains the plot number as well as the time point that a litter bag",
        "from the plot was collected. There are 2 time point columns, called",
        "'samp' and 'timePoint'. The values correspond so that 1 -> 0, 2 -> 3",
        "3 -> 5, 4 -> 6, where the value to the left of each arrow is from",
        "'samp' while the value to the right of the arrow is from",
        "'timePoint'. The column 'compositionParameter' describes the",
        "specific parameter that describes microbial community composition.",
        "The parameters called 'bac' and 'fun', I believe, correspond to",
        "bacterial and fungal reads, respectively. 'fb' is calculated by",
        "dividing fungal reads by bacterial reads and is used to make the",
        "figure on fungal:bacterial ratios. 'shan' is Shannon's taxonomic",
        "diversity. The column 'value' contains the actual data value of",
        "each composition parameter for each plot and time point."]
readMe = pd.DataFrame({"Note": note})
exportPath = folderPath/"Community composition.xlsx"
if os.path.exists(exportPath) is False:
    with pd.ExcelWriter(exportPath) as w:
        readMe.to_excel(w, "ReadMe", index=False)
        combined.to_excel(w, "Data", index=False)
