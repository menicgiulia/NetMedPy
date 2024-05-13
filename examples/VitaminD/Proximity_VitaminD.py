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

    query = "select distinct diseaseNID,geneNID,associationType,sentence,score,EL,EI from geneDiseaseNetwork where score >= 0.3"
    query = "select * from (" + query + ") as gd_net natural join "
    query = query + " (select geneNID,geneName,geneDescription from geneAttributes) natural join "
    query = query + " (select * from diseaseAttributes)"

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

    if disease=='None':
        return (None,len(chol_targets),np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    else:
        disease_genes = gda.query("`Disease name` == @disease")
        disease_genes = set(disease_genes['Gene Symbol'])

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


if __name__=='__main__':
    ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    with open(ppi_file, 'rb') as file:
        ppi = pickle.load(file)

    ppi = netmedpy.extract_lcc(list(ppi.nodes), ppi)
    ppi_nodes = set(ppi.nodes)

    ##########################Para Wang Genes##############################
    # gda = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/associations/disgenet_wang.xlsx"
    # gda = pd.read_excel(gda)
    # gda = gda[gda['Gene Symbol'].isin(set(ppi_nodes))]
    #################################################################3


    #######################Para DISGENET#####################################
    # dg_file = "D:/data/disgenet/disgenet_2020.db"
    # gda = load_DisGeNet(dg_file)
    # gda.rename(columns={'geneName':'Gene Symbol','diseaseName':'Disease name'},inplace=True)
    # gda = gda[gda['Gene Symbol'].isin(set(ppi_nodes))]
    ##########################################################################

    # disease = gda[gda['Disease name'].str.contains("factor vii",case=False)]


    #########################Para deisy GDA########################################
    dg_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/associations/GDA_Filtered_04042022.csv"

    gda = pd.read_csv(dg_file)
    gda = gda.query("Strong > 0 | Weak> 0")
    gda = gda[gda['HGNC_Symbol'].isin(ppi_nodes)]
    gda = gda.rename(columns={'HGNC_Symbol': 'Gene Symbol', 'NewName':'Disease name'})
    ###############################################################################

    disease = gda[gda['Disease name'].str.contains("Triple Negative Breast Neoplasms",case=False)]

    ################Diseases Disgenet and wang#############################
    # diseases = ['Inflammation',
    #             'Coronary Artery Disease',
    #             'Coronary heart disease',
    #             'Lung Neoplasms',
    #             # 'Breast Neoplasms',
    #             'Triple Negative Breast Neoplasms',
    #             'Brain Neoplasms',
    #             'Pulmonary Emphysema',
    #             'alpha-Thalassemia',
    #             'Essential Hypertension',
    #             'Pancreatic Neoplasm',
    #             'Rett Syndrome, Atypical',
    #             'Rubinstein-Taybi Syndrome',
    #             'Pulmonary Hypertension',
    #             'Prader-Willi Syndrome',
    #             'Carcinoma',
    #             'Asthma',
    #             'Cerebrovascular accident',
    #             'Factor VII Deficiency',
    #             'beta Thalassemia',
    #             'Huntington Disease',
    #             'Fragile X Syndrome',
    #             'Colorectal Neoplasms',
    #             'Cystic Fibrosis',
    #             'Muscular Dystrophy, Duchenne',
    #             'None'
    #             ]

    #########################Diseases Channing###############################
    # diseases = ["inflammation",
    #          "coronary artery disease",
    #          "coronary heart disease",
    #          "lung neoplasms",
    #          "breast neoplasms",
    #          "brain neoplasms",
    #          "COPD",
    #          "alpha-Thalassemia", #
    #          "hypertension",
    #          "pancreatic neoplasms",
    #          "Rett_syndrome", #
    #          "RubinsteinTaybi_syndrome", #
    #          "pulmonary hypertension",
    #          "PraderWilli_syndrome", #
    #          "carcinoma",
    #          "asthma",
    #          "cerebrovascular disease",
    #          "FactorVII_deficiency", #
    #          "beta thalassemia",
    #          "huntington disease",
    #          "FragileX_syndrome", #
    #          "colorectal neoplasms",
    #          "cystic fibrosis",
    #          "muscular dystrophy duchenne"
    #          "FactorIX_deficiency" #
    #          ]

    #####################################



    ########################Diseases Deisy######################################
    diseases = ['inflammation',
                'coronary artery disease',
                'coronary disease',
                'lung neoplasms',
                # 'Breast Neoplasms',
                'breast neoplasms',
                'brain neoplasms',
                'pulmonary emphysema',
                'None',
                'hypertension',
                'pancreatic neoplasms',
                'rett syndrome',
                'None',
                'hypertension pulmonary',
                'prader willi syndrome',
                'carcinoma',
                'asthma',
                'None',
                'None',
                'beta thalassemia',
                'huntington disease',
                'None',
                'colorectal neoplasms',
                'cystic fibrosis',
                'muscular dystrophy duchenne',
                'None'
                ]



    dGenes = number_disease_genes(diseases, gda)




    ######################Channing drug targets#################################
    targets_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"
    with open(targets_file,'rb') as file:
        targets= pickle.load(file)

    d,t = targets
    chol_targets = set(t['Cholecalciferol'])
    chol_targets = chol_targets & ppi_nodes


    ####################Drugbank#######################################
    # drugbank = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/Drug_Bank/DB_Drug_Targets_2023.csv"
    # drugbank = pd.read_csv(drugbank)
    # drugbank = drugbank[drugbank['Gene_Target'].isin(ppi_nodes)]

    # chol_targets = drugbank.query("Name=='Cholecalciferol'")
    # chol_targets = set(chol_targets.Gene_Target)


    print("Calculating distance matrix")
    # dMatrix = netmedpy.all_pair_distances(ppi, distance='shortest_path',n_processors=19,n_tastks=100)
    # netmedpy.save_distances(dMatrix, "D:/data/matDist/channing/distance_unweighted_channing.pkl")
    dMatrix = netmedpy.load_distances("D:/data/matDist/channing/distance_unweighted_channing.pkl")


    print("Evaluating Screening")
    results = screening(diseases,ppi,gda,chol_targets,dMatrix)

    results.to_csv("C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/screening_deisy_channingDT.csv",index=False)



    #######################################Plot results#######################################
    def plot_results():
        main_dir = "data/vitamin_d/"

        df = pd.read_csv(main_dir + "screening_channing.csv")
        df_z = df[['Disease','Z-score']]
        df_z.columns = ['Disease','Channing']

        df_dg = df[['Disease','Disease genes']]
        df_dg.columns = ['Disease','Channing']

        df= pd.read_csv(main_dir + "screening_disgenet_channingDT.csv")
        df_z['DisGeNet_chDT'] = df['Z-score']
        df_dg['DisGeNet_chDT'] = df['Disease genes']
        df_z['Disease'] = df['Disease']
        df_dg['Disease'] = df['Disease']

        df= pd.read_csv(main_dir + "screening_deisy_channingDT.csv")
        df_z['Deisy_chDT'] = df['Z-score']
        df_dg['Deisy_chDT'] = df['Disease genes']

        df= pd.read_csv(main_dir + "screening_disgenet_drugbankDT.csv")
        df_z['DisGeNet_DB'] = df['Z-score']
        df_dg['DisGeNet_DB'] = df['Disease genes']

        df= pd.read_csv(main_dir + "screening_deisy_drugbankDT.csv")
        df_z['Deisy_DB'] = df['Z-score']
        df_dg['Deisy_DB'] = df['Disease genes']
        
        

        # Plotting
        plt.figure(figsize=(10, 7))
        pal = plt.rcParams['axes.prop_cycle'].by_key()['color']

        df = df_z.iloc[:-1,]
        plt.axhline(y=0,color='gray',linewidth=2,linestyle='--')
        # Plotting both 'ch' and 'dg' as line plots with points
        plt.plot(range(len(df['Disease'])), df['DisGeNet_DB'], marker='o', 
                  color=pal[3], linestyle='--', label='DisGeNet_DB')
        plt.plot(range(len(df['Disease'])), df['Deisy_DB'], marker='o', 
                  color=pal[4],linestyle='--', label='Deisy_DB')

        
        plt.plot(range(len(df['Disease'])),df['DisGeNet_chDT'], 
                  marker='o', color=pal[2],label='DisGeNet chDT')
        plt.plot(range(len(df['Disease'])), df['Deisy_chDT'], 
                  marker='o', color=pal[1], label='Deisy chDT')
        plt.plot(range(len(df['Disease'])),df['Channing'], marker='o', 
                 linewidth=2,color='black',label='Channing')

        # Rotating and aligning the x-axis labels
        plt.xticks(ticks = range(len(df['Disease'])),
                   labels = df['Disease'],rotation=45, ha='right',fontsize=12)

        plt.ylabel('Z-Score')
        plt.legend(loc='lower right')

        plt.tight_layout()
        plt.show()
