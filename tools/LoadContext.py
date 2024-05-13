# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 15:43:56 2023

@author: aalda
"""

import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

plt.style.use('seaborn-v0_8-colorblind')
plt.rcParams["figure.figsize"] = [9, 6]
plt.rcParams["figure.autolayout"] = True

plt.rcParams['font.size'] = 15
plt.rcParams.update({'font.size':15})
# Set the axes labels font size
plt.rc('axes', labelsize=15)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=15)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=15)



class Context:

    def __init__(self):
        PPI_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/networks/PPI/PPI_2022_v2.csv"
        PPI_df = pd.read_csv(PPI_file)

        self.PPI = nx.from_pandas_edgelist(PPI_df,source="HGNC_Symbol.1",target="HGNC_Symbol.2")


        GDA = pd.read_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/associations/GDA_Filtered_04042022.csv")
        self.GDA = GDA.query('Weak >= 1 or Strong >= 1 or Incompatible >= 1')


        self.DISTANCE_FILE = "C:/Users/aalda/OneDrive/Documentos/ppi_all_distances/all_distances.pkl"
