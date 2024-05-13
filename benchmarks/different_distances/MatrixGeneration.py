# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 18:12:19 2023

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



if __name__=="__main__":
    c= Cronometer.Cronometer()

    ###############Random Networks###########################

    #m0 = 5  # Initial number of nodes
    #m = 2   # Number of edges to attach from a new node to existing nodes

    # Generate the graph
    #G_ba = nx.barabasi_albert_graph(18000, m, seed=42)




    #Rename nodes
    # nodes = list(G_ba.nodes())

    # mapping = { n:f"n{n}" for n in nodes}

    # G = nx.relabel_nodes(G_ba,mapping)
    ###################################################

    ####################Real network##############################
    context = LoadContext.Context()

    G = context.PPI
    #############################################################

    metrics = ["shortest_path","random_walk","biased_random_walk","communicability"]

    output_path = "D:/data/ppi_distance_matrices/"

    times = []

    for m in metrics:
        c.tick()
        D = network_metrics.all_pair_distances(G, distance=m,reset = 0.2)
        c.tock()

        print(f"Finished {m}")
        print(f"Time {c.format_seconds()}")


        outfile = f"{output_path}distances_{m}.pkl"
        network_metrics.save_distances(D, outfile)

        times.append(c.elapsed_seconds)

        del D

    df = pd.DataFrame(columns=metrics)
    df.loc[0] = times

    df.to_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/all_pair_distances/ppi_matrix_times.csv",index=False)
