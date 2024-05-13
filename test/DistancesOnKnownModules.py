# -*- coding: utf-8 -*-
"""
Created on Wed Oct 4 18:44:39 2023
@author: Andres Aldana Gonzalez (a.aldana@northeastern.edu)

This test creates 4 scale free modules, connects them in sequence and evaluates
proximity and separation between them.
"""

import networkx as nx
import netmedpy.NetMedPy as metrics
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


if __name__=="__main__":

    G = nx.Graph()


    edges = [ ('A','B'),
             ('A','C'),
             ('B','C'),
             ('B','E'),
             ('C','D'),
             ('D','E'),
             ('D','F'),
             ('F','G'),
             ('E','G'),
             ('E','H'),
             ('G','I'),
             ('H','I'),
             ('H','J'),
             ('I','J'),
             ('J','K'),
             ('I','L'),
             ('K','L'),
            ]

    G.add_edges_from(edges)

    plt.figure(figsize=(9,6))
    nx.draw(G,
            with_labels=True,
            node_color='skyblue',
            node_size=1000,
            edge_color='gray')
    plt.show()


    mods = [{'A','B','C','D','E'},
            {'D','E','F','G'},
            {'H','I','J'},
            {'J','K','L'},
            {'F','G','H','I'},
           ]
    n_modules = len(mods)




    # --- CALCULATING DISTANCES ---

    # Calculate the all-pair shortest path distances in G
    D = metrics.all_pair_distances(G)

    # Save and then reload the distance matrix
    metrics.save_distances(D, "distances.pkl")
    D = metrics.load_distances("distances.pkl")



    #--- PROXIMITY ANALYSIS ---

    #Calculate the proximity between pairs of modules

    prox_mat_z = []
    prox_mat_p = []
    prox_mat_raw = []

    for i in range(n_modules):
        g1 = mods[i]
        l_z = []
        l_p = []
        l_r = []
        for j in range(n_modules):
            print(f"\rProximity {i},{j}                   ",end="")
            g2 = mods[j]
            p = metrics.proximity(G,g1,g2,D,n_iter=3000,degree_preserving='exact',bin_size=100)
            l_z.append(p['z_score'])
            l_p.append(p['p_val'])
            l_r.append(p['raw_amspl'])

        prox_mat_z.append(l_z)
        prox_mat_p.append(l_p)
        prox_mat_raw.append(l_r)

    df_z = pd.DataFrame(prox_mat_z)
    df_p = pd.DataFrame(prox_mat_p)
    df_raw = pd.DataFrame(prox_mat_raw)

    # # --- SEPARATION ANALYSIS ---

    sep_mat_z = []
    sep_mat_p = []
    sep_mat_raw = []

    for i in range(n_modules):
        g1 = mods[i]
        l_z = []
        l_p = []
        l_r = []
        for j in range(n_modules):
            print(f"\rSeparation {i},{j}                   ",end="")
            g2 = mods[j]
            p = metrics.separation_z_score(G,g1,g2,D,n_iter=3000,degree_preserving='exact',bin_size=100)
            l_z.append(p['z_score'])
            l_p.append(p['p_val'])
            l_r.append(p['raw_separation'])

        sep_mat_z.append(l_z)
        sep_mat_p.append(l_p)
        sep_mat_raw.append(l_r)

    df_sep_z = pd.DataFrame(sep_mat_z)
    df_sep_p = pd.DataFrame(sep_mat_p)
    df_sep_raw = pd.DataFrame(sep_mat_raw)
