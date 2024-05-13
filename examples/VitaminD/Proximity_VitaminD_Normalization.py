

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

import seaborn as sns
import tools.plot_styles


@ray.remote
def process_disease(disease, ppi, gda, chol_targets, distance_matrix):
    dlist,dgenes = gda
    disease_genes = set(dgenes[disease])

    prox = netmedpy.proximity(ppi, chol_targets, disease_genes, distance_matrix,n_iter=5000)

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

def plot_group(norm):
    norm = "N2"
    base_path = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/normalizationTest/"

    ####################Test N1#####################
    files = ["shortest_path",f"random_walk_{norm}", f"biased_random_walk_{norm}",f"communicability_{norm}"]

    df = pd.DataFrame()

    for f in files:
        if df.empty:
            df = pd.read_csv(f"{base_path}{f}.csv")
            #df = df[["Disease","Z-score"]]
            df = df[["Disease","amspl"]]
            df.columns = ["Disease",f"{f}"]
        else:
            df2 = pd.read_csv(f"{base_path}{f}.csv")
            #df2 = df2[["Disease","Z-score"]]
            df2 = df2[["Disease","amspl"]]
            df2.columns = ["Disease",f"{f}"]

            df = pd.merge(df,df2,on="Disease")

    df = df.sort_values(by='shortest_path',ascending=True)


    plt.figure(figsize=(9,6))
    pal = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.axhline(y=0,color='gray',linewidth=2,linestyle='--')


    plt.xticks(ticks = range(len(df['Disease'])),
               labels = df['Disease'],rotation=45, ha='right',fontsize=12)

    plt.plot(range(len(df['Disease'])), df[f'random_walk_{norm}'], marker='o',
              color=pal[0], linestyle='-', linewidth = 1, label=f'Random Walk {norm}')

    plt.plot(range(len(df['Disease'])), df[f'biased_random_walk_{norm}'], marker='o',
              color=pal[1], linestyle='-', linewidth = 1, label=f'Biased Random Walk {norm}')

    plt.plot(range(len(df['Disease'])), df[f'communicability_{norm}'], marker='o',
              color=pal[2], linestyle='-', linewidth = 1, label=f'Communicability {norm}')

    plt.plot(range(len(df['Disease'])), df['shortest_path'], marker='o',
              color="black", linestyle='-', linewidth = 2, label='Shortest Path')

    plt.ylabel("Proximity Z-Score")
    plt.legend(loc='best',frameon=False)
    plt.tight_layout()
    plt.grid(visible=None)
    plt.show()


    df_cor = df.drop("Disease",axis=1)
    corr_matrix = df_cor.corr()

    plt.figure(figsize=(9,6))
    sns.set(font_scale=1.5)
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', linewidths=.5)
    plt.show()


    #######################Rank plots
    cols = df.columns
    cols = cols[1:]


    df_rank = pd.DataFrame()
    for c in cols:

        dfNew = df[["Disease",c]]
        dfNew = dfNew.sort_values(by = c)
        rank = list(range(len(dfNew.index)))
        dfNew["rank"] = rank

        dfNew = dfNew[["Disease","rank"]]
        dfNew.columns = ["Disease",c]

        if df_rank.empty:
            df_rank = dfNew
        else:
            df_rank = pd.merge(df_rank,dfNew,on="Disease")


    plt.figure(figsize=(12,7))
    plt.xticks(ticks = range(len(df['Disease'])),
               labels = df_rank['Disease'],rotation=45, ha='right',fontsize=12)
    plt_index = 0
    for c in cols:
        plt.plot(range(len(df_rank['Disease'])),df_rank[c], label= c, marker='o',
              color=pal[plt_index], linestyle='-', linewidth = 1)
        plt_index = plt_index + 1
    plt.legend(loc="best")
    plt.ylabel("Rank")
    plt.grid(visible=None)
    plt.ylim(0,30)
    plt.show()



    df_cor = df_rank.drop("Disease",axis=1)
    corr_matrix = df_cor.corr(method="spearman")

    plt.figure(figsize=(9,6))
    sns.set(font_scale=1.5)
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', linewidths=.5)
    plt.show()




def plot_N1N2(metric):
    base_path = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/normalizationTest/"

    ####################Test N1#####################
    files = ["shortest_path",f"{metric}_N1", f"{metric}_N2"]

    df = pd.DataFrame()

    for f in files:
        if df.empty:
            df = pd.read_csv(f"{base_path}{f}.csv")
            df = df[["Disease","Z-score"]]
            df.columns = ["Disease",f"{f}"]
        else:
            df2 = pd.read_csv(f"{base_path}{f}.csv")
            df2 = df2[["Disease","Z-score"]]
            df2.columns = ["Disease",f"{f}"]

            df = pd.merge(df,df2,on="Disease")

    df = df.sort_values(by='shortest_path',ascending=True)


    plt.figure(figsize=(9,6))
    pal = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.axhline(y=0,color='gray',linewidth=2,linestyle='--')


    plt.xticks(ticks = range(len(df['Disease'])),
               labels = df['Disease'],rotation=45, ha='right',fontsize=12)

    plt.plot(range(len(df['Disease'])), df[f'{metric}_N1'], marker='o',
              color=pal[0], linestyle='-', linewidth = 1, label=f'N1')

    plt.plot(range(len(df['Disease'])), df[f'{metric}_N2'], marker='o',
              color=pal[1], linestyle='-', linewidth = 1, label=f'N2')

    plt.plot(range(len(df['Disease'])), df['shortest_path'], marker='o',
              color="black", linestyle='-', linewidth = 2, label='Shortest Path')

    plt.ylabel("Proximity Z-Score")
    plt.legend(loc='best',frameon=False)
    plt.title(f"{metric}")
    plt.tight_layout()
    plt.show()





def plot_results():
    plot_group("N1")
    plot_group("N2")

    plot_N1N2("random_walk")
    plot_N1N2("biased_random_walk")
    plot_N1N2("communicability")



def communicability_null_models(ppi, distance_matrix, chol_targets, disease_genes):
    control_file = f"C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/normalizationTest/random_walk_N2.csv"
    control = pd.read_csv(control_file)
    control = control.sort_values(by="Z-score")

    df = control[["Disease","Z-score"]]
    df.columns = ["Disease","Control"]

    null_models = ["degree_match","log_binning","uniform"]
    ctargets = {"Cholecalciferol":chol_targets}

    for n in null_models:
        print(f"Evaluating {n}")
        sc = netmedpy.screening(ctargets,disease_genes,ppi,distance_matrix,score="proximity",prop="z_score",
                                null_model=n,n_iter=10000,bin_size=100,symmetric=False,n_procs=None)

        sc = sc.T
        sc = sc.reset_index()
        sc = sc.rename(columns={'index':'Disease','Cholecalciferol':n})

        df = pd.merge(df,sc,on="Disease")

    df = df.sort_values(by="Control")

    return df


def plot_communicability_null_models(df):


    plt.figure()
    pal = plt.rcParams['axes.prop_cycle'].by_key()['color']
    pal = ["black"] + pal

    plt.axhline(y=0,color='gray',linewidth=2,linestyle='--')


    plt.xticks(ticks = range(len(df['Disease'])),
               labels = df['Disease'],rotation=45, ha='right',fontsize=12)


    plot_cols = list(df.columns)
    plot_cols.remove('Disease')

    for i,c in enumerate(plot_cols):
        plt.plot(range(len(df['Disease'])), df[c], marker='o',
                  color=pal[i], linestyle='-', linewidth = 1, label=c)


    plt.legend(loc='best',frameon="False")
    plt.ylabel("Proximity Z-Score")
    plt.show()
    
    
    df = df_null_models
    df_cor = df.drop("Disease",axis=1)
    
    for c in df_cor.columns:
        df_cor[c] = pd.to_numeric(df_cor[c])
        
    corr_matrix = df_cor.corr(method="spearman")
    

    plt.figure(figsize=(9,6))
    sns.set(font_scale=1.5)
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', linewidths=.5)
    plt.show()
    
    
    


if __name__=='__main__':

    matrix_file= "communicability_N2"

    ppi_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    with open(ppi_file, 'rb') as file:
        ppi = pickle.load(file)

    ppi = netmedpy.extract_lcc(list(ppi.nodes), ppi)
    ppi_nodes = set(ppi.nodes)


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


    targets_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/chem_targets_10192023.pkl"
    with open(targets_file,'rb') as file:
        targets= pickle.load(file)

    d,t = targets
    chol_targets = set(t['Cholecalciferol'])
    chol_targets = chol_targets & ppi_nodes

    print("Calculating distance matrix")
    dMatrix = netmedpy.load_distances(f"D:/data/matDist/channing/ppi932023/{matrix_file}.pkl")


    print("Evaluating Screening")
    #results = screening(dlist,ppi,gene_set,chol_targets,dMatrix)

    #results.to_csv(f"C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/NetMedPy/data/vitamin_d/normalizationTest/{matrix_file}.csv",index=False)


    #plot_results()

    df_null_models = communicability_null_models(ppi,dMatrix,chol_targets,dgenes)
    plot_communicability_null_models(df_null_models)
    
    
    
