# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:04:43 2022

@author: dalensis (minor modifications from Brian Chung)
This script is copied from GitHub user dalensis on this pingouin issue:
https://github.com/raphaelvallat/pingouin/issues/205

Minor modifications are made to make this look cleaner and also remove certain
unnecessary lines, such as an if statement that checks if the type of the
p-value is a number or a string. The purpose of this script is to construct a
"compact letter display", essentially letters that are assigned to treatment
groups in a multiple comparisons test to designate groups that are similar or
different to each other.
"""
import string
import pandas as pd
import os
from pathlib import Path

#launch the main function with the output of the pairwise comparison and the CI (for example 99 for alpha=0.01)


def main(df, CI):
    if len(df.index)<2: df = df.rename(columns = {"p-unc" : "pval"})    #the pval column  has different names based on test and numerosity
    else: df = df.rename(columns = {"p-corr" : "pval"})
    
    groups = sorted(set(df["A"].unique()).union((set(df["B"].unique())))) #take all the names from columns A and B 
    letters = list(string.ascii_uppercase)[:len(groups)]
    cldgroups = letters
    
    #the following algoritm is a semplification of the classical cld, 
    
    cld = pd.DataFrame(list(zip(groups, letters, cldgroups)))
    for row in df.itertuples():                                 
        if not type(df["pval"][row[0]]) is str and df["pval"][row[0]]>(1-CI/100):
            cld.iat[groups.index(df["A"][row[0]]), 2] +=  cld.iat[groups.index(df["B"][row[0]]), 1]
            cld.iat[groups.index(df["B"][row[0]]), 2] +=  cld.iat[groups.index(df["A"][row[0]]), 1]
    
    cld[2] = cld[2].apply(lambda x: "".join(sorted(x)))
    
   #this part will reassign the final name to the group, for sure there are more elegant way of doing it
    cld = cld.sort_values(cld.columns[2], key=lambda x: x.str.len())
    cld["groups"]=""
    letters = list(string.ascii_uppercase)
    unique = []
    for item in cld[2]:

        for fitem in cld["groups"].unique():
            for c in range (0, len(fitem)):
                if not set(unique).issuperset(set(fitem[c])):
                    unique.append(fitem[c])
        g=len(unique)
        
        for kitem in cld[1]:
            if kitem in item:
                if cld["groups"].loc[cld[1]==kitem].iloc[0]=="": cld["groups"].loc[cld[1]==kitem]+=letters[g]
                if not len(set(cld["groups"].loc[cld[1]==kitem].iloc[0]).intersection(cld.loc[cld[2]==item,"groups"].iloc[0]))>0:
                    if letters[g] not in list(cld["groups"].loc[cld[1]==kitem].iloc[0]): cld["groups"].loc[cld[1]==kitem]+=letters[g]
                    if letters[g] not in list(cld["groups"].loc[cld[2]==item].iloc[0]): cld["groups"].loc[cld[2]==item]+=letters[g]
   
    cld = cld.sort_values("groups") 
    print(cld)
    return(cld) #return the df. In my base script i catch it, save to xls, and use the groups to tag the column of the plot.


# Importing alkane, timePoint x Vegetation raw Tukey results with columns that
# have been modified in Excel
cwd = Path(os.getcwd())
filesNfolders = os.listdir(cwd)
testDataPath = cwd/'alkane, timePoint x Vegetation.xlsx'
testData = pd.read_excel(testDataPath)
testData.pval = testData["pval"].astype(float)
testOutput = main(testData, 95)

'''Ok looks like this algorithm works. Fuck yes.'''
