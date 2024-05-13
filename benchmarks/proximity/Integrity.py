# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 14:55:22 2023

@author: aalda
"""
import proximity.distances as old_distances
import dnetmedpy.NetMedPy as new_distances
import networkx as nx
import numpy as np
import tools.LoadContext as LoadContext
import tools.Cronometer as Cronometer
import pandas as pd
import matplotlib.pyplot as plt
import random
from proximity.network import Network as agNetwork
import tools.General as genTools
import multiprocessing
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

import gseapy as gp


global_context = None

def evaluate_disease_pair(CONTEXT):
    diseases = genTools.get_disease_names(CONTEXT)

    d1 = random.choice(diseases)
    d2 = random.choice(diseases)

    dg1 = genTools.get_disease_genes(d1,CONTEXT)
    dg2 = genTools.get_disease_genes(d2, CONTEXT)

    ppi = CONTEXT.PPI

    try:
        graph = agNetwork(ppi)
        pr = old_distances.proximity(graph, dg1, dg2,n_iter=1000)
        pr_1 = pr['z_score']
    except Exception as e:
        pr_1 = np.nan


    try:
        pr = new_distances.proximity_allCalc(ppi, dg1, dg2,'exact',10,0)
        pr_2 = pr['z_score']
    except Exception as e:
        pr_2 = np.nan

    if pr_1 != np.nan and pr_2 != np.nan:
        dif = pr_1-pr_2
    else:
        dif = np.nan

    return [d1,d2,pr_1,pr_2, dif]



def all_pair_distances(graph):
    global graph_global
    graph_global = graph

    # Create a Pool of workers
    # The initializer and initargs are used to set the graph_global for each process
    with Pool(processes=cpu_count(), initializer=init_pool, initargs=(graph,)) as pool:
        # Use a generator instead of a list for the nodes
        results = pool.map(worker, (node for node in graph.nodes()))

    return {node: distances for node, distances in results}





def init_pool(CONTEXT):
    global global_context
    global_context = CONTEXT


def worker(k):
    res = evaluate_disease_pair(global_context)
    return res

# def eval_k_pairs(k, CONTEXT, outFile):
#     global global_context
#     global_context = CONTEXT

#     pbar = tqdm(total=k*2)


#     with Pool(processes = cpu_count(), initializer=init_pool, initargs=(CONTEXT,)) as pool:
#         results = pool.map(worker,list(range(k)))

#     pbar.close()

#     df = pd.DataFrame(data=results,columns=['D1','D2','z1','z2','z1-z2'])


#     df.to_csv(outFile,index=False)
#     return results




def eval_k_pairs(k,CONTEXT,outFile):
   global global_context
   global_context = CONTEXT

   with tqdm(total=k, desc="Processing", dynamic_ncols=True, position=0, leave=True) as pbar:
       with Pool(processes=cpu_count(), initializer=init_pool, initargs=(CONTEXT,)) as pool:
           results = []
           for result in pool.imap_unordered(worker, list(range(k))):
               results.append(result)
               pbar.update()

   df = pd.DataFrame(data=results, columns=['D1', 'D2', 'z1', 'z2', 'z1-z2'])
   return df


def plot_series(file):

    df = pd.read_csv(file)

    x = list(range(len(df.index)))

    plt.figure()
    plt.scatter(x, df.z1, marker='o',s=40,label='Aproximate DPR',color="blue")
    plt.scatter(x, df.z2, marker='o',s=40,label='Exact DPR',color="red")
    plt.xlabel("Disease pairs")
    plt.ylabel("Proximity")
    plt.legend(loc = 'best',frameon=False,markerscale=1.0)
    plt.show()

    plt.figure()
    plt.hist(df['z1-z2'], bins=50)
    plt.xlabel('z_aprox - z_exact')
    plt.ylabel("Count")
    plt.show()


def plot_z_vs_jaccard(file,CONTEXT):
    df = pd.read_csv(file)

    corz1 = df['z1'].corr(df['jaccard'])
    corz2 = df['z2'].corr(df['jaccard'])


    plt.figure()
    plt.scatter(df.z1, df.jaccard, marker='o',s=40,label=f"Aproximate DPR (P={format(corz1,'.2f')})",color="blue")
    plt.scatter(df.z2, df.jaccard, marker='o',s=40,label=f"Exact DPR (P={format(corz2,'.2f')})",color="red")
    plt.xlabel("Proximity")
    plt.ylabel("Jaccard")
    plt.legend(loc = 'best',frameon=False,markerscale=1.0)
    # plt.title(f"Pz1 = {format(corz1,'.2f')} Pz2={format(corz2,'.2f')}")
    plt.show()



def common_pathways(name1,name2, gene_set ,CONTEXT):
    bg = list(CONTEXT.PPI.nodes)

    gs = 'GO_Biological_Process_2023' if gene_set == 'BP' else 'Reactome_2022'


    d1genes = genTools.get_disease_genes(name1, CONTEXT)
    d2genes = genTools.get_disease_genes(name2, CONTEXT)

    bps1 = gp.enrichr(gene_list=list(d1genes),gene_sets=gs
                      ,organism='Human',background=bg)
    bps2 = gp.enrichr(gene_list=list(d2genes),gene_sets=gs
                      ,organism='Human',background=bg)


    bps1 = bps1.results.query(" `Adjusted P-value`< 0.05 ")
    bps2 = bps2.results.query(" `Adjusted P-value`< 0.05 ")


    terms1 = set(bps1.Term)
    terms2 = set(bps2.Term)

    return len(terms1 & terms2)/len(terms1 | terms2)


def complement_bp(df,gene_set,CONTEXT):
    jaccard = []

    colName = 'bp_'+gene_set

    for i in df.index:
        print(f"Processing {i+1}/{len(df.index)}")
        r = df.iloc[i]

        d1 = r['D1']
        d2 = r['D2']

        c = common_pathways(d1, d2, gene_set, CONTEXT)

        jaccard.append(c)

    df[colName] = jaccard



def plot_z_bp(file,gene_set,yAxis,CONTEXT):

    df = pd.read_csv(file)

    colName = 'bp_'+gene_set
    corz1 = df['z1'].corr(df[colName])
    corz2 = df['z2'].corr(df[colName])


    plt.figure()
    plt.scatter(df.z1, df[colName], marker='o',s=40,label=f"Aproximate DPR (P={format(corz1,'.2f')})",color="blue")
    plt.scatter(df.z2, df[colName], marker='o',s=40,label=f"Exact DPR (P={format(corz2,'.2f')})",color="red")
    plt.xlabel("Proximity")
    plt.ylabel(yAxis)
    plt.legend(loc = 'best',frameon=False,markerscale=1.0)
    # plt.title(f"Pz1 = {format(corz1,'.2f')} Pz2={format(corz2,'.2f')}")
    plt.show()



def complement_with_jaccards(file,CONTEXT):
    df = pd.read_csv(file)

    df = df.dropna()
    df.reset_index(drop=True, inplace=True)

    df['abs'] = df['z1-z2'].abs()

    ###########Pure jaccard
    jaccard = []

    for i in df.index:

        r = df.iloc[i]

        d1 = r['D1']
        d2 = r['D2']


        d1_genes = genTools.get_disease_genes(d1,CONTEXT)
        d2_genes = genTools.get_disease_genes(d2,CONTEXT)

        inters = set(d1_genes) & set(d2_genes)
        un = set(d1_genes) | set(d2_genes)

        j = len(inters)/ len(un)

        jaccard.append(j)

    df['jaccard'] = jaccard


    ###########Biological Processes

    complement_bp(df, 'BP', CONTEXT)

    #############Reactome Pathways
    complement_bp(df, 'PW', CONTEXT)

    df.to_csv(file,index=False)



if __name__=='__main__':

    file= "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/proximity/integrity/integrity.csv"

    CONTEXT = LoadContext.Context()

    k = 100

    # print("Initializing pool")

    # c = Cronometer.Cronometer()

    # c.tick()
    # results = eval_k_pairs(k, CONTEXT,file)
    # c.tock()

    # print("Elapsed time: " + c.format_seconds())

    # print("Complementing dataframe with jaccard distances")
    # complement_with_jaccards(file, CONTEXT)
    # print("Done")

    df = plot_series(file)

    plot_z_vs_jaccard(file, CONTEXT)

    plot_z_bp(file, 'BP', 'Jaccard - GO Processes', CONTEXT)

    plot_z_bp(file, 'PW', 'Jaccard - Reactome Pathways', CONTEXT)
