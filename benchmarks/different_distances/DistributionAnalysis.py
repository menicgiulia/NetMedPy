# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 23:33:11 2023

@author: aalda
"""



import networkx as nx
import tools.Cronometer as Cronometer
import numpy as np
import netmedpy.DistanceMatrix as dMatrix
import netmedpy.NetMedPy as network_metrics
import matplotlib.pyplot as plt
import tools.LoadContext as LoadContext
import pandas as pd
import seaborn as sns


if __name__=="__main__":
    c= Cronometer.Cronometer()

    context = LoadContext.Context()

    #############First we plot the time to generate the matrices##
    df = pd.read_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/all_pair_distances/ppi_matrix_times.csv")

    times = df.loc[0,:]
    labels = ["Shortest Path","RWR","Biased RWR","Communicability"]

    plt.figure()
    plt.bar(labels,times/60,width=0.6)
    plt.ylabel("Time (m)")
    plt.suptitle("Matrix generation time")
    plt.show()

    ######################Matrix Distribution plots
    matrix_path = "D:/data/ppi_distance_matrices/"

    D = network_metrics.load_distances(f"{matrix_path}distances_communicability.pkl")

    values = D.matrix.flatten()
    sample = np.random.choice(values, size=10000, replace=False)

    max(values)

    sample = [v for v in sample if v < 0.2]

    plt.figure()
    plt.hist(sample, bins=100,align='mid')
    # plt.suptitle('Value Distribution')
    plt.xlabel('Communicability')
    plt.ylabel('Frequency')
    plt.yscale("log")
    plt.show()
