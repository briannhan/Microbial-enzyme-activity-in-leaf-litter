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


def makeCLD(df, alpha=0.05):
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
    df.dropna(axis="columns")
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
    cld[3] = ""

    for row in df.itertuples():
        if df["p-adj"][row[0]] > (alpha):
            cld.iat[groups.index(df["group1"][row[0]]), 2] += cld.iat[groups.index(df["group2"][row[0]]), 1]
            cld.iat[groups.index(df["group2"][row[0]]), 2] += cld.iat[groups.index(df["group1"][row[0]]), 1]

        if df["p-adj"][row[0]] < (alpha):
            cld.iat[groups.index(df["group1"][row[0]]), 3] += cld.iat[groups.index(df["group2"][row[0]]), 1]
            cld.iat[groups.index(df["group2"][row[0]]), 3] += cld.iat[groups.index(df["group1"][row[0]]), 1]

    cld[2] = cld[2].apply(lambda x: "".join(sorted(x)))
    cld[3] = cld[3].apply(lambda x: "".join(sorted(x)))
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

                # Checking if there are forbidden pairing (proposition of solution to the imperfect script)
                if kitem in ' '.join(cld[3][cld["labels"] == letters[g]]):
                    g = len(unique)+1

                # Checking if columns 1 & 2 of cld share at least 1 letter
                if len(set(cld["labels"].loc[cld[1] == kitem].iloc[0]).intersection(cld.loc[cld[2] == item, "labels"].iloc[0])) <= 0:
                    if letters[g] not in list(cld["labels"].loc[cld[1] == kitem].iloc[0]):
                        cld["labels"].loc[cld[1] == kitem] += letters[g]
                    if letters[g] not in list(cld["labels"].loc[cld[2] == item].iloc[0]):
                        cld["labels"].loc[cld[2] == item] += letters[g]

    cld = cld.sort_values("labels")
    print(cld)
    print('\n')
    cld.drop(columns=[1, 2, 3], inplace=True)
    print(cld)
    print('\n')
    print('\n')
    return(cld)


# Importing alkane, timePoint x Vegetation raw Tukey results with columns that
# have been modified in Excel. This is to test the function
# cwd = Path(os.getcwd())
# filesNfolders = os.listdir(cwd)
# testDataPath = cwd/'Vmax, timePoint x Vegetation.xlsx'
# testData = pd.read_excel(testDataPath)
# testOutput = makeCLD(testData, 0.05)

'''Ok looks like this algorithm works. Fuck yes.

Nope. It does not work on NAG Vmax time x veg. The raw results states that
5 x CSS and 6 x CSS are different, but this algorithm makes it so that they
both share the letter 'b', so it's not perfect. Goddamn it.

Ok so the guy made some changes to it. Let's see if it works now.

Tested it again on NAG Vmax time x veg, and it works this time. Fuck yes.
'''


def checkCLD(data, cld):
    """
    Tests whether the compact letter display matches the raw Tukey results.

    Parameters
    ----------
    data : Pandas dataframe
        Contains the raw Tukey test results.
    cld : Pandas dataframe
        Contains the compact letter display that is created from the raw Tukey
        test results.

    Returns
    -------
    A string, the values of which are either 'Good' or 'Bad' to denote whether
    the compact letter display is correct or wrong when compared to the Tukey
    test results.

    """
    # Splitting the letter strings into individual characters
    cldCopy = cld.copy()
    for index, row in cldCopy.iterrows():
        letterStr = row.labels
        cldCopy.loc[index, "labels"] = list(letterStr)  # Splitting string by
        # converting to list, which creates a list of individual characters

    # Checking to see if the compact letter display matches the raw Tukey
    # results
    errorStr = ""
    for index, row in data.iterrows():
        testResult = row.reject
        group1 = row.group1
        group2 = row.group2
        group1cldRow = cldCopy[cldCopy["groups"] == group1]
        group1cld = group1cldRow["labels"].tolist()[0]
        group1cld = set(group1cld)
        group2cldRow = cldCopy[cldCopy["groups"] == group2]
        group2cld = group2cldRow["labels"].tolist()[0]
        group2cld = set(group2cld)
        intersection = group1cld & group2cld  # Set operation that creates an
        # intersection of the 2 sets

        # Test the scenario where the 2 groups are similar to each other, with
        # the "reject" column from the raw Tukey results being "FALSE". In this
        # scenario, they should share at least 1 letter. If they don't, then
        # print an error message.
        if testResult == "FALSE":
            if len(intersection) == 0:
                errorStr = "The cld does not match raw Tukey test results"
                resultStr = "Bad"
                print(errorStr)
                break

        # Test the scenario where the 2 groups are different from each other,
        # with the "reject" column from the raw Tukey results being "TRUE".
        # They should not share any letters. If they do share at least 1
        # letter, then print an error message
        elif testResult == "TRUE":
            if len(intersection) > 0:
                errorStr = "The cld does not match raw Tukey test results"
                resultStr = "Bad"
                print(errorStr)
                break
    if errorStr != "The cld does not match raw Tukey test results":
        print("Compact letter display matches raw Tukey results")
        resultStr = "Good"

    return resultStr


# check = checkCLD(testData, testOutput)
# myManualCLD = pd.read_excel("Vmax, timePoint x Vegetation, correct.xlsx")
# myManualCLD.dropna(axis=1, inplace=True)
# check = checkCLD(testData, myManualCLD)


def main(testResults, alpha=0.05):
    """
    Creates and checks the compact letter display by running the 2 functions
    above. The intended use of this function/method is that when the user
    imports this module into another script, the user simply needs to invoke
    this function, which would then output the compact letter display from raw
    Tukey results. The user doesn't need to call any of the above 2 methods to
    create the compact letter display

    Parameters
    ----------
    testResults : Pandas dataframe
        The raw Tukey test results as formatted by statsmodels.
    alpha : float, optional
        The alpha value that the user used in their Tukey test. The default is
        0.05.

    Returns
    -------
    A Pandas dataframe that represents the compact letter display.

    """
    cld = makeCLD(testResults)
    cldTest = checkCLD(testResults, cld)
    if cldTest == "Good":
        return cld
    elif cldTest == "Bad":
        print(cldTest)
    return


# Importing alkane, timePoint x Vegetation raw Tukey results with columns that
# have been modified in Excel. This is to test the function
# cwd = Path(os.getcwd())
# filesNfolders = os.listdir(cwd)
# testDataPath = cwd/'Vmax, timePoint x Vegetation.xlsx'
# testData = pd.read_excel(testDataPath)
# testOutput = main(testData, 0.05)
# testOutput = makeCLD(testData, 0.05)
