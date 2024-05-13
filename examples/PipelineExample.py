# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:34:10 2024

@author: aalda
"""
import netmedpy.DistanceMatrix as dmat
import netmedpy.NetMedPy as metrics

import pandas as pd
import networkx as nx
import pickle

if __name__=='__main__':

    #Loading PPI
    file_path = "C:/Users/aalda/CCNR Dropbox/Biology/Aldana Gonzalez, Andres/Gordana/Interaction Networks/ppi932023.pkl"

    with open(file_path,'rb') as file:
        ppi = pickle.load(file)

    ppi = metrics.extract_lcc(ppi.nodes, ppi)

    #Disease genes
    file_gda = "C:/Users/aalda/CCNR Dropbox/Biology/Aldana Gonzalez, Andres/Gordana/Gene Disease Associations/disgenet_wang.xlsx"
    gda = pd.read_excel(file_gda)


    #Load DrugBank
    drugbank_file = "C:/Users/aalda/CCNR Dropbox/Biology/Aldana Gonzalez, Andres/Gordana/DrugBank/Processed_data/all_drugbank_drugs.csv"
    drugbank = pd.read_csv(drugbank_file)
    drugbank = drugbank.query("Status=='approved' & organism=='Humans'")


    #Drugs related to asthma
    drugbank_asthma = drugbank[drugbank["Indication"].str.contains("asthma",case=False,na=False)]
    drugbank_asthma = drugbank_asthma[drugbank_asthma["Gene_Target"].isin(set(ppi.nodes))]


    #drugs related to Alzheimer
    drugbank_ah = drugbank[drugbank["Indication"].str.contains("alzheimer",case=False,na=False)]
    drugbank_ah = drugbank_ah[drugbank_ah["Gene_Target"].isin(set(ppi.nodes))]


    #Select drug targets from one drug related to asthma
    fluti = drugbank_asthma.query("Name=='Fluticasone propionate'")
    fluti_targets = set(fluti.Gene_Target)

    #Select drug targets from one drug related to Alzheimer
    galan = drugbank_ah.query("Name=='Galantamine'")
    galan_targets = set(galan.Gene_Target)


    #Pick Asthma disease genes
    asthma_dg = gda[gda["Disease name"]=="Asthma"]
    asthma_dg = set(asthma_dg["Gene Symbol"])

    asthma_dg = asthma_dg.intersection(set(ppi.nodes))


    lcc = metrics.extract_lcc(asthma_dg, ppi)
    lcc_nodes = list(lcc.nodes)
    stat = metrics.lcc_significance(ppi, lcc_nodes,n_iter=2000)


    D = metrics.all_pair_distances(ppi, distance="shortest_path")

    #metrics.save_distances(D, "Distances.pkl")
    #D = metrics.load_distances("Distances.pkl")


    p_fluti = metrics.proximity(ppi, fluti_targets, asthma_dg, D,null_model='degree_match',n_iter=2000,symmetric=True)


    p_galan = metrics.proximity(ppi, galan_targets, asthma_dg, D,null_model='degree_match',n_iter=2000,symmetric=True)


    p_fluti['z_score']
    p_galan['z_score']


    p_fluti['p_value_single_tail']
    p_galan['p_val_single_tail']

    p_fluti['raw_amspl']
    p_galan['raw_amspl']
