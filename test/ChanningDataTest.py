# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 20:42:40 2024

@author: aalda
"""

import pickle as pkl
import pandas as pd
import numpy as np
import networkx as nx
import netmedpy.DistanceMatrix as DistanceMatrix
import netmedpy.NetMedPy as metrics
import matplotlib.pyplot as plt
import ray
from multiprocessing import cpu_count
import seaborn as sns


@ray.remote
def sep_row(i,disList,Distance, geneDict,ppi):

    row = []
    d1Genes = set( geneDict[ disList[i] ] )

    for j in range(len(disList)):
        if i == j:
            row.append(0)
        else:
            print(f"Sep({disList[i]},{disList[j]}) ")
            d2Genes = set( geneDict[ disList[j] ] )
            s = metrics.separation_z_score(ppi, d1Genes, d2Genes, Distance,
                                           null_model="degree_match",
                                           n_iter=1000)

            row.append(s)

    #print(f"Finished {i}/{len(disList)} ")

    return (i,row)


def generate_sep_results(disList,Distance,geneDict,ppi):
    ray.shutdown()
    ray.init(num_cpus=cpu_count())

    disList_ref = ray.put(disList)
    Distance_ref = ray.put(Distance)
    geneDict_ref = ray.put(geneDict)
    ppi_ref = ray.put(ppi)


    results = [sep_row.remote(i,disList_ref,Distance_ref,geneDict_ref,ppi_ref)
            for i in range(len(disList))]

    results = ray.get(results)

    ray.shutdown()


    return results


def transform_to_DataFrame(results,disList):
    mat = np.zeros((len(disList),len(disList)))

    for (i,row) in results:

        ar = []

        for j in range(len(row)):
            elem = row[j]
            if i == j:
                ar.append(elem)
            else:
                ar.append(elem['z_score'])

        mat[i] = np.array(ar)

    df = pd.DataFrame(mat)

    df.columns = disList
    df.index = disList

    return df



def plot_heatmap(matrix,):
    plt.figure(figsize=(15,9))  # Adjust the size as needed

    # Create a heatmap
    ax = sns.heatmap(matrix, annot=True, fmt=".2f", cmap="bwr", center=0, linewidths=.5, cbar_kws={"shrink": .75})
    ax.tick_params(axis='both', which='major', labelsize=13)
    plt.show()



if __name__=="__main__":

    f_Mppi = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/PPI202207.txt"
    f_Jppi = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    f_targets = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"
    f_diseases = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_10192023.pkl"


    Mppi = pd.read_csv(f_Mppi,sep ='\t',header=None,
                       dtype={0:np.int32,1:str,2:np.int32,3:str})
    Mppi.columns = ["Gene A","Symbol A","Gene B","Symbol B"]
    Mppi =Mppi[["Symbol A","Symbol B"]]
    Mppi = nx.from_pandas_edgelist(Mppi,source="Symbol A",target="Symbol B")



    with open(f_Jppi, 'rb') as file:
        Jppi = pd.read_pickle(file)
    Jppi = metrics.extract_lcc(Jppi.nodes, Jppi)
    ppi_nodes = set(Jppi.nodes())

    with open(f_targets, 'rb') as file:
        targets = pd.read_pickle(file)
    medlist,drug_targets = targets



    with open(f_diseases, 'rb') as file:
        diseases = pd.read_pickle(file)
    dlist, dgenes = diseases

    for d in dlist:
        l =  dgenes[d]
        l = set(l) & ppi_nodes

        dgenes[d] = l
    diseases = (dlist,dgenes)


    Jppi_lcc = metrics.extract_lcc(list(Jppi.nodes), Jppi)
    print(f"LCC={len(Jppi_lcc)} Total={len(Jppi)}")

    Jppi = Jppi_lcc

    distFile = "D:/data/matDist/channing/distance_unweighted.pkl"

    D = metrics.all_pair_distances(Jppi, distance="shortest_path")
    #metrics.save_distances(D, distFile)

    # D = metrics.load_distances(distFile)

    res = generate_sep_results(diseases[0], D, diseases[1], Jppi)

    df = transform_to_DataFrame(res, diseases[0])

    plot_heatmap(df)

row = res[1][1][3]['z_score']
