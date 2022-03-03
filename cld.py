# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:04:43 2022

@author: dalensis (minor modifications from Brian Chung)
This module is copied from GitHub user dalensis on this pingouin issue:
https://github.com/raphaelvallat/pingouin/issues/205

Minor modifications are made to make this look cleaner and also remove certain
unnecessary lines, such as an if statement that checks if the type of the
p-value is a number or a string. The purpose of this module is to construct a
"compact letter display", essentially letters that are assigned to treatment
groups in a multiple comparisons test to designate groups that are similar or
different to each other.
"""
import string
import pandas as pd
import os
from pathlib import Path


def main(df, alpha=0.05):
    '''
    Creates a compact letter display. This creates a dataframe consisting of
    2 columns, a column containing the treatment groups and a column containing
    the letters that have been assigned to the treatment groups. These letters
    are part of what's called the compact letter display. Treatment groups that
    share at least 1 letter are similar to each other, while treatment groups
    that don't share any letters are significantly different from each other.

    Parameters
    ----------
    df : Pandas dataframe
        Pandas dataframe containing raw Tukey test results from statsmodels.
    alpha : float
        The alpha value. The default is 0.05.

    Returns
    -------
    A dataframe representing the compact letter display, created from the Tukey
    test results.

    '''
    df["p-adj"] = df["p-adj"].astype(float)

    # Creating a list of the different treatment groups from Tukey's
    group1 = set(df.group1.tolist())  # Dropping duplicates by creating a set
    group2 = set(df.group2.tolist())  # Dropping duplicates by creating a set
    groupSet = group1 | group2  # Set operation that creates a union of 2 sets
    groups = sorted(list(groupSet))

    # Creating lists of letters that will be assigned to treatment groups
    letters = list(string.ascii_lowercase)[:len(groups)]
    cldgroups = letters

    # the following algoritm is a simplification of the classical cld,

    cld = pd.DataFrame(list(zip(groups, letters, cldgroups)))
    for row in df.itertuples():
        if df["p-adj"][row[0]] > (alpha):
            cld.iat[groups.index(df["group1"][row[0]]), 2] += cld.iat[groups.index(df["group2"][row[0]]), 1]
            cld.iat[groups.index(df["group2"][row[0]]), 2] += cld.iat[groups.index(df["group1"][row[0]]), 1]

    cld[2] = cld[2].apply(lambda x: "".join(sorted(x)))
    cld.rename(columns={0: "groups"}, inplace=True)

    # this part will reassign the final name to the group
    # for sure there are more elegant ways of doing this
    cld = cld.sort_values(cld.columns[2], key=lambda x: x.str.len())
    cld["labels"] = ""
    letters = list(string.ascii_lowercase)
    unique = []
    for item in cld[2]:

        for fitem in cld["labels"].unique():
            for c in range(0, len(fitem)):
                if not set(unique).issuperset(set(fitem[c])):
                    unique.append(fitem[c])
        g = len(unique)

        for kitem in cld[1]:
            if kitem in item:
                if cld["labels"].loc[cld[1] == kitem].iloc[0] == "":
                    cld["labels"].loc[cld[1] == kitem] += letters[g]

                # Checking if columns 1 & 2 of cld share at least 1 letter
                if len(set(cld["labels"].loc[cld[1] == kitem].iloc[0]).intersection(cld.loc[cld[2] == item, "labels"].iloc[0])) <= 0:
                    if letters[g] not in list(cld["labels"].loc[cld[1] == kitem].iloc[0]):
                        cld["labels"].loc[cld[1] == kitem] += letters[g]
                    if letters[g] not in list(cld["labels"].loc[cld[2] == item].iloc[0]):
                        cld["labels"].loc[cld[2] == item] += letters[g]

    cld = cld.sort_values("labels")
    print(cld)
    print('\n')
    cld.drop(columns=[1, 2], inplace=True)
    print(cld)
    print('\n')
    print('\n')
    return(cld)


# Importing alkane, timePoint x Vegetation raw Tukey results with columns that
# have been modified in Excel. This is to test the function
cwd = Path(os.getcwd())
filesNfolders = os.listdir(cwd)
testDataPath = cwd/'alkane, timePoint x Vegetation.xlsx'
testData = pd.read_excel(testDataPath)
testOutput = main(testData, 0.05)

'''Ok looks like this algorithm works. Fuck yes.

Nope. It does not work on NAG Vmax time x veg. The raw results states that
5 x CSS and 6 x CSS are different, but this algorithm makes it so that they
both share the letter 'b', so it's not perfect. Goddamn it. '''
