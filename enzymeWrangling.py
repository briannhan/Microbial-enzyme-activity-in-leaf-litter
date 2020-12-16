# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 18:51:39 2020

@author: Brian Nhan Thien Chung
This is a module that is meant to facilitate the wrangling of enzyme activity
data and to help plot them.
"""
import pandas as pd

# %%
# Functions to wrangle both hydrolytic & oxidative enzyme data


def processLongNamesAndWells(enzymeData):
    """Processes the long sample names & wells in an enzyme dataframe that
    occurred from reading in a dataframe of raw enzyme data. The long sample
    names are in a column, which are then split up into 2 columns. The Well
    column will also be used to create 2 new columns


    Parameters
    ----------
    enzymeData : Pandas dataframe
        Dataframe of raw enzyme data from reading in a raw enzyme data text
        file. Has 3 columns: Well, Long sample name, Assay. Well will be used
        to create 2 new columns, while Long sample name will also be used
        to create 2 new columns

    Returns
    -------
    The dataframe with the names & wells processed & with new columns

    """
