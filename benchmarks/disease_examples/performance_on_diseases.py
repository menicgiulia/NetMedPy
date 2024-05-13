# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:35:32 2023

@author: aalda
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tools.LoadContext as LoadContext
import netmedpy.DistanceMatrix as DistanceMatrix
import netmedpy.NetMedPy as NetworkMetrics
import networkx as nx
import gseapy as gp
import seaborn as sns
import ray



def lcc(graph, nodes):
    subgraph = graph.subgraph(nodes)
    connected_components = list(nx.connected_components(subgraph))
    largest_component = max(connected_components, key=len)
    return subgraph.subgraph(largest_component)


def jaccard(A,B):
    a = set(A)
    b = set(B)

    return len(a & b)/len(a | b)


if __name__=="__main__":


    context = LoadContext.Context()

    ppi = context.PPI
    gda = context.GDA


    genes = gda.query("NewName=='alzheimer disease'").HGNC_Symbol
    alz = lcc(ppi,genes)
    alz = list(alz.nodes)

    genes = gda.query("NewName=='arthritis rheumatoid'").HGNC_Symbol
    art = lcc(ppi,genes)
    art = list(art.nodes)

    genes = gda.query("NewName=='crohn disease'").HGNC_Symbol
    crohn = lcc(ppi,genes)
    crohn = list(crohn.nodes)


    ##################Enrichment Analysis#####################
    dis_list = [alz,art,crohn]
    paths = []

    for l in dis_list:
        enr = gp.enrichr(gene_list=l,gene_sets='Reactome_2022',organism='Human')
        # enr = gp.enrichr(gene_list=l,gene_sets='GO_Biological_Process_2023',organism='Human')
        enr = enr.results
        enr = enr.query("`Adjusted P-value` < 0.05")
        paths.append(list(enr.Term))


    jacs = []
    for i in range(len(paths)):
        ne = (i+1)%len(paths)

        j = jaccard(paths[i], paths[ne])
        jacs.append(j)

    labels = ["AH-RA","RA-CH","AH-CH"]


    plt.figure()
    plt.bar(labels,jacs,width=0.6)
    plt.ylabel("Jaccard Index",fontsize=16)
    plt.suptitle("Overlap between common Reactome pathways")
    plt.show()


    ##############Calculate separation and produce matrix
    path = "D:/data/ppi_distance_matrices/"
    types = ["shortest_path","random_walk","biased_random_walk","communicability"]
    labels = ["Shortest Path","RWR","Biased RWR","Communicability"]

    pairs = [ (alz,art),(art,crohn), (alz,crohn)]
    pairs_labels = ["AH-RA","RA-CH","AH-CH"]

    distances = np.zeros((len(types),len(pairs)))


    for i in range(len(types)):
        file = f"{path}distances_{types[i]}.pkl"
        D = NetworkMetrics.load_distances(file)

        for j in range(len(pairs)):
            a,b = pairs[j]

            print(f"Calculating {pairs_labels[j]} for {labels[i]} ")
            s = NetworkMetrics.separation_z_score(ppi, a, b,D,
                                                  degree_preserving="exact",n_iter=5000)

            distances[i,j] = s['z_score']



    # # Plotting the bar chart
    palette = sns.color_palette(palette='YlOrBr')
    palette = palette[(len(palette)-len(pairs_labels)):]

    fig, ax = plt.subplots()
    bar_width = 0.2
    opacity = 1.0

    for i, pair in enumerate(pairs_labels):
        ax.bar(np.arange(len(labels)) + bar_width * i, distances[:, i], bar_width,
                alpha=opacity, label=pair,color = palette[i])

    # ax.set_xlabel('Distance Definitions')
    ax.set_ylabel('Separation',fontsize=18)
    ax.set_title('Separation between pairs of diseases')
    ax.set_xticks(np.arange(len(labels)) + bar_width)
    ax.set_xticklabels(labels,fontsize=18)
    ax.legend(loc="best")

    plt.show()
