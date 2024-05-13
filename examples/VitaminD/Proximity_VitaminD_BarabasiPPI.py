# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 16:28:57 2024

@author: aalda
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 18:49:34 2024

@author: aalda
"""
import networkx as nx
import numpy as np
import pandas as pd

import netmedpy.NetMedPy as netmedpy
import matplotlib.pyplot as plt

import pickle

import ray

import sqlite3 as sq



def load_DisGeNet(file):
    conn = sq.connect(file)
    cur = conn.cursor()

    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")

    tables = cur.fetchall()

    query = 'select * from (select * from geneDiseaseNetwork where score >= 0.3) as gd_net natural join'
    query = query + ' (select geneNID,geneId,geneName,geneDescription from geneAttributes) natural join'
    query = query + ' diseaseAttributes'

    diseaseGenes = pd.read_sql(query,conn)

    return diseaseGenes


def number_disease_genes(diseases, gda):
    dg = []
    dg_contains = []

    for d in diseases:
        q = gda[gda['Disease name'].str.contains(d,case=False)]
        dg_contains.append(len(q.index))

        q = gda.query("`Disease name`==@d")
        dg.append(len(set(q['Gene Symbol'])))

    return pd.DataFrame({'Disease':diseases,'Equal_DG':dg,'Contains_DG':dg_contains})


@ray.remote
def process_disease(disease, ppi, gda, chol_targets, distance_matrix):

    if disease =='None':
        return (None,len(chol_targets),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    else:
        dlist,dgenes = gda
        disease_genes = set(dgenes[disease])

        prox = netmedpy.proximity(ppi, chol_targets, disease_genes, distance_matrix,n_iter=10000)

        print(f"{disease} finished")

        return (disease, len(chol_targets), len(disease_genes),prox['z_score'], prox['raw_amspl'],
                prox['p_value_single_tail'], prox['p_value_double_tail'], prox['d_mu'],prox['d_sigma'])


def screening(disease_list,ppi,gda,chol_targets,distance_matrix):
    ray.shutdown()
    ray.init()

    ppi_ref = ray.put(ppi)
    gda_ref = ray.put(gda)
    chol_targets_ref = ray.put(chol_targets)
    distance_matrix_ref = ray.put(distance_matrix)

    res = [process_disease.remote(disease, ppi_ref, gda_ref,chol_targets_ref,distance_matrix_ref)
           for disease in disease_list]

    results = ray.get(res)

    ray.shutdown()

    columns = ['Disease','Drug Targets','Disease genes','Z-score','amspl','p_value_st','p_value_dt','mu','sigma']

    res = pd.DataFrame(results,columns = columns)

    return res


def plot_DB_vs_CPE():
    base_dir = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/barabasiPPI/"

    uCPE = pd.read_csv(base_dir + "screening_unweighted_chanDT.csv")
    uCPE = uCPE[['Disease','Z-score']]
    uCPE.columns = ['Disease','z_cpe']

    uDB = pd.read_csv(base_dir + "screening_unweighted_DrugBank.csv")
    uDB = uDB[['Disease','Z-score']]
    uDB.columns = ['Disease','z_db']

    wCPE = pd.read_csv(base_dir + "screening_weighted_cpe.csv")
    wCPE = wCPE[['Disease','Z-score']]
    wCPE.columns = ['Disease','z_wcpe']

    df = pd.merge(uCPE,uDB,how='inner',on='Disease')
    df = pd.merge(df,wCPE,how='inner',on='Disease')

    df = df.sort_values(by='z_cpe')


    # Plotting
    plt.figure(figsize=(10, 7))
    pal = plt.rcParams['axes.prop_cycle'].by_key()['color']

    plt.axhline(y=0,color='gray',linewidth=2,linestyle='--')
    # Plotting both 'ch' and 'dg' as line plots with points
    plt.plot(range(len(df['Disease'])), df['z_cpe'], marker='o',
              color="black", linestyle='-', label='CPE')
    plt.plot(range(len(df['Disease'])), df['z_db'], marker='o',
              color=pal[0],linestyle='-', label='Drug Bank')
    plt.plot(range(len(df['Disease'])), df['z_wcpe'], marker='o',
              color=pal[1],linestyle='-', label='Weighted CPE')


    # Rotating and aligning the x-axis labels
    plt.xticks(ticks = range(len(df['Disease'])),
               labels = df['Disease'],rotation=45, ha='right',fontsize=12)

    plt.ylabel('Z-Score')
    plt.legend(loc='best')

    plt.tight_layout()
    plt.show()



if __name__=='__main__':
    # ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    # with open(ppi_file, 'rb') as file:
    #     ppi = pickle.load(file)

    # ppi = netmedpy.extract_lcc(list(ppi.nodes), ppi)


    # ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/networks/PPI/PPI_2023-05-10.csv"
    # ppi = pd.read_csv(ppi_file)
    # ppi = nx.from_pandas_edgelist(ppi,source='HGNC_Symbol.1',target='HGNC_Symbol.2')

    filename = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/PPIimprovement_and_RA_RNAseq/data/weightedNets/PPI_2023-05-10_kendall_GSE97810.pkl"
    with open(filename, 'rb') as file:
        ppi = pickle.load(file)


    ppi_nodes = set(ppi.nodes)


    channin_disease_genes_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/disease_set_01262024.pkl"
    with open(channin_disease_genes_file, 'rb') as file:
        gene_set = pickle.load(file)
    dlist, dgenes = gene_set

    for d in dlist:
        l =  dgenes[d]
        l = set(l) & ppi_nodes

        dgenes[d] = l


    diseases = ["inflammation",
            "coronary artery disease",
            "coronary heart disease",
            "lung neoplasms",
            "breast neoplasms",
            "brain neoplasms",
            "COPD",
            "alpha-Thalassemia", #
            "hypertension",
            "pancreatic neoplasms",
            "Rett_syndrome", #
            "RubinsteinTaybi_syndrome", #
            "pulmonary hypertension",
            "PraderWilli_syndrome", #
            "carcinoma",
            "asthma",
            "cerebrovascular disease",
            "FactorVII_deficiency", #
            "beta thalassemia",
            "huntington disease",
            "FragileX_syndrome", #
            "colorectal neoplasms",
            "cystic fibrosis",
            "muscular dystrophy duchenne"
            ]


    #d_registered = set(dlist) & set(diseases)
    gene_set = (dlist,dgenes)


    ######################FOR CHANNING DRUG TARGETS#########################
    targets_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"
    with open(targets_file,'rb') as file:
        targets= pickle.load(file)

    d,t = targets
    chol_targets = set(t['Cholecalciferol'])
    chol_targets = chol_targets & ppi_nodes


    #####################3 FOR DRUGBANK ###########################3
    # drugbank = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/Drug_Bank/DB_Drug_Targets_2023.csv"
    # drugbank = pd.read_csv(drugbank)
    # drugbank = drugbank[drugbank['Gene_Target'].isin(ppi_nodes)]
    # drugbank = drugbank.query("organism=='Humans'")

    # chol_targets = drugbank.query("Name=='Cholecalciferol'")
    # chol_targets = set(chol_targets.Gene_Target)


    print("Calculating distance matrix")
    dMatrix = netmedpy.all_pair_distances(ppi, distance='shortest_path',n_processors=19,n_tastks=100)
    netmedpy.save_distances(dMatrix, "D:/data/matDist/barabasiPPI/distance_weighted_kendall_GSE_PPI_2023-05-10.pkl")
    #dMatrix = netmedpy.load_distances("D:/data/matDist/barabasiPPI/distance_unweighted_PPI_2023-05-10.pkl")


    print("Evaluating Screening")
    results = screening(diseases,ppi,gene_set,chol_targets,dMatrix)
    results.to_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/barabasiPPI/screening_weighted_cpe.csv",index=False)
