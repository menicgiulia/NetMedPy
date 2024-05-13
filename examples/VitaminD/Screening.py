# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:43:00 2024

@author: aalda
"""

import networkx as nx
import numpy as np
import pandas as pd

import netmedpy.NetMedPy as netmedpy
import matplotlib.pyplot as plt

import pickle
import ray

import seaborn as sns

def merge_lung_neoplasm():
    channin_disease_genes_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_01262024.pkl"
    with open(channin_disease_genes_file, 'rb') as file:
        gene_set = pickle.load(file)
    dlist, dgenes = gene_set

    lneo = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/lungneoplasms.csv"
    lneo = pd.read_csv(lneo)

    lneogenes = set(lneo['lungneoplasms'])

    dgenes['lung neoplasms'] = lneogenes

    outFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_03142024.pkl"

    out = (dlist,dgenes)

    with open(outFile, 'wb') as file:
        pickle.dump(out,file)


def create_drug_dictionary():
    chol_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"

    with open(chol_file,'rb') as file:
        d = pickle.load(file)

    lista,dictionary = d

    co = dictionary['Cholecalciferol']

    dt_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/vitaminD_drugs/cpe.csv"

    dt = pd.read_csv(dt_file)

    num_pchembl = pd.to_numeric(dt['ave_pchembl'], errors='coerce').notna()

    dt_filter = dt[num_pchembl].copy()
    dt_filter['ave_pchembl'] = pd.to_numeric(dt_filter['ave_pchembl'])
    dt_filter = dt_filter.query("ave_pchembl> 3")

    dt = dt_filter

    names = dt.Prescription.unique()


    finalDict = {'Cholecalciferol':co}

    for n in names:
        drug_targets = dt.query("Prescription==@n")
        drug_targets = set(drug_targets.hgnc_symbol)

        finalDict[n] = drug_targets

    return finalDict

# @ray.remote
# def process_disease_targets(disease, drug, ppi, diseases, drugs, distance_matrix):
#     drug_targets = drugs[drug]
#     disease_genes = diseases[disease]

#     prox = netmedpy.proximity(ppi, drug_targets, disease_genes, distance_matrix,n_iter=10000)

#     print(f"{disease}-{drug} finished")

#     return (disease, drug,prox['z_score'], prox['raw_amspl'],
#             prox['p_value_single_tail'], prox['p_value_double_tail'], prox['d_mu'],prox['d_sigma'])



# def screening(disease_list,disease_genes, drug_targets, ppi, distance_matrix):
#     ray.shutdown()
#     ray.init()

#     ppi_ref = ray.put(ppi)
#     dl_ref = ray.put(disease_genes)
#     dt_ref = ray.put(drug_targets)
#     distance_matrix_ref = ray.put(distance_matrix)

#     futures = []

#     for disease in disease_genes:
#         for drug in drug_targets:
#             future = process_disease_targets.remote(disease,drug,
#                                                     ppi_ref, dl_ref,dt_ref,
#                                                     distance_matrix_ref)
#             futures.append(future)

#     results = ray.get(futures)

#     ray.shutdown()


#     columns = ["Cholecalciferol","Ibuprofen","Dexamethasone", "Atorvastatin",
#                "Amlodipine","Salmeterol","Fluticasone furoate","Pralsetinib","Nintedanib",
#                "Risperidone","Fluoxetine"]


#     df = pd.DataFrame(columns=columns, index=disease_list)

#     for r in results:
#         dis = r[0]
#         drug = r[1]
#         prox = r[2]

#         df.loc[dis,drug] = prox

#     df = df.sort_values(by='Cholecalciferol',ascending=True)

#     return df


def plot_matrix():
    df = pd.read_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/screening_all_drugs.csv",index_col=0)

    plt.figure(figsize=(9, 5))  # Adjust the figure size as needed
    sns.set_context("notebook", font_scale=1.25)

    sns.heatmap(df, annot=False, fmt="g", cmap='seismic',center=0,vmin=-10,vmax=2)

    # Rotate and align the x-axis labels
    plt.xticks(rotation=45, ha='right',fontsize=16)
    plt.yticks(fontsize=16)

    # Show the plot
    plt.show()

if __name__=='__main__':

    ## Load PPI
    ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    with open(ppi_file, 'rb') as file:
        ppi = pickle.load(file)

    ppi = netmedpy.extract_lcc(list(ppi.nodes), ppi)
    ppi_nodes = set(ppi.nodes)


    ## Load disease genes
    channin_disease_genes_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_03142024.pkl"
    with open(channin_disease_genes_file, 'rb') as file:
        gene_set = pickle.load(file)
    dlist, dgenes = gene_set

    for d in dlist:
        l =  dgenes[d]
        l = set(l) & ppi_nodes

        dgenes[d] = l

    disease_list = ["inflammation","coronary artery disease","COPD",
                "lung neoplasms","FragileX_syndrome"]

    disease_genes = {}
    for d in disease_list:
        disease_genes[d] = dgenes[d]


    ## Load drug_targets
    drug_targets = create_drug_dictionary()

    for k,v in drug_targets.items():
        v = set(v) & ppi_nodes
        drug_targets[k] = v


    print("Calculating distance matrix")
    #dMatrix = netmedpy.all_pair_distances(ppi, distance='shortest_path',n_processors=19,n_tastks=100)
    #netmedpy.save_distances(dMatrix, "D:/data/matDist/channing/ppi932023/shortest_path.pkl")
    dMatrix = netmedpy.load_distances("D:/data/matDist/channing/ppi932023/shortest_path.pkl")


    print("Evaluating Screening")
    #results = screening(disease_list, disease_genes, drug_targets,ppi,dMatrix)
    results = netmedpy.screening(drug_targets, disease_genes, ppi, dMatrix,score="proximity",prop="z_score",
                                 n_iter=100,n_procs=19)


    #results.to_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/screening_all_drugs.csv",index=True)
