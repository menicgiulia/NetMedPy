# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 16:22:23 2024

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
import tools.plot_styles








if __name__=='__main__':

    data_path = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/examples/VitaminD/paper/data/"

    ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    with open(ppi_file, 'rb') as file:
        ppi = pickle.load(file)

    ppi = netmedpy.extract_lcc(list(ppi.nodes), ppi)
    ppi_nodes = set(ppi.nodes)


    with open(data_path + "ppi_network.pkl","wb") as file:
        pickle.dump(ppi,file)


    channin_disease_genes_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_03142024.pkl"
    with open(channin_disease_genes_file, 'rb') as file:
        gene_set = pickle.load(file)
    dlist, dgenes = gene_set

    delete_list = ["AutoCore","AutoCore_Interleukin and cytokine","AutoCore_NFkB and chromatin","AutoCore_Proteasome complex",
                   "AutoCore_Transcription factors, DNA binding complex","dnmts literature","epigenetic modifiers",
                   "epigenetic modifiers DNA","epigenetic modifiers DNA meth","one carbon pathway"]

    for d in delete_list:
        del dgenes[d]
    dlist = dgenes.keys()

    for d in dlist:
        l =  dgenes[d]
        l = set(l) & ppi_nodes

        dgenes[d] = l

    gene_set = (dlist,dgenes)


    alias = pd.read_csv(data_path + "Alias.csv",index_col=0)
    r_dgenes = {}

    for k,v in dgenes.items():
        r_dgenes[alias.loc[k,"correct_name"]] = v


    with open(data_path + "disease_genes.pkl",'wb') as file:
        pickle.dump(r_dgenes,file)



    targets_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"
    with open(targets_file,'rb') as file:
        targets= pickle.load(file)

    d,t = targets
    chol_targets = set(t['Cholecalciferol'])
    chol_targets = chol_targets & ppi_nodes

    with open(data_path + "vitd_targets.pkl", 'wb') as file:
        pickle.dump(chol_targets,file)
