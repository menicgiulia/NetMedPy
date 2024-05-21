# -*- coding: utf-8 -*-
"""
Exploring Vitamin D’s Impact on Autoimmune, Cardiovascular, and Cancer Diseases: A Network Medicine Perspective Case Study withNetMedPy.

@author: Andres Aldana Gonzalez
@Date: 04-25-2024

This script evaluates the role of Vitamin D in the modulation of
autoimmune diseases, cardiovascular diseases and cancer from a
network medicine perspective and reproduces the results presented in the paper "Title of the paper"

To run this script install packages networkx and ray with

pip install networkx
pip install ray
"""

#Set the working directory to VitaminD directory
#import os
#os.chdir("/user_path_to/NetMedPy/examples/VitaminD/")


import networkx as nx
import numpy as np
import pandas as pd

import netmedpy.NetMedPy as netmedpy
import matplotlib.pyplot as plt

import pickle
import ray

import seaborn as sns


plt.style.use('seaborn-v0_8-colorblind')
plt.rcParams["figure.figsize"] = [11,7]
plt.rcParams["figure.autolayout"] = True
font = 20

plt.rcParams['font.size'] = font
plt.rcParams.update({'font.size':font})
# Set the axes labels font size
plt.rc('axes', labelsize=font)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=font)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=font)


def plot_lcc_significance(lcc_data):

    df = lcc_data.sort_values(by="size",ascending=False)
    df = df.fillna(0)



    #plot size
    plt.figure()
    plt.bar(df['disease'], df['size'])
    plt.axhline(y=10,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC Size')
    plt.xticks(rotation=45, ha='right')  # Rotate labels 45 degrees and align them right
    plt.yscale("log")
    plt.show()

    #plot z-score
    plt.figure()
    plt.bar(df['disease'], df['zscore'])
    plt.axhline(y=2,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC Z-Score')
    plt.xticks(rotation=45, ha='right')  # Rotate labels 45 degrees and align them right
    plt.show()

    #plot p-value
    plt.figure()
    plt.bar(df['disease'], df['pval'])
    plt.axhline(y=0.05,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC p-value')
    plt.xticks(rotation=45, ha='right')  # Rotate labels 45 degrees and align them right
    plt.yscale("log")
    plt.show()


def plot_histograms(inflammation, factorix):

    plt.figure()
    sns.kdeplot(inflammation['dist'], color='blue', fill=True,
                alpha=0.5, label='Inflammation')
    sns.kdeplot(factorix['dist'], color='red', fill=True,
                alpha=0.5, label='Factor IX Def.')


    plt.axvline(x=inflammation['raw_amspl'],linewidth=10,
                color='blue',alpha=0.5)
    plt.text(x=inflammation['raw_amspl'], y=8,
             s=f"  Z = {inflammation['z_score']:.2f}",
             color='blue', verticalalignment='top',
             horizontalalignment='left',fontsize=20)


    plt.axvline(x=factorix['raw_amspl'],linewidth=10,color='red',alpha=0.5)

    plt.text(x=factorix['raw_amspl'], y=8,
             s=f"  Z = {factorix['z_score']:.2f}",
             color='red', verticalalignment='top',
             horizontalalignment='left',fontsize=20)

    plt.xlabel('AMSPL from Vitamin D')
    plt.ylabel('Density')
    plt.ylim(0,9)
    plt.legend(loc='upper center',frameon=False)
    plt.show()



def plot_screening(sdata):
    vdata = sdata.T
    vdata = vdata.reset_index()
    vdata.columns = ["disease","zscore"]
    vdata = vdata.sort_values(by="zscore",ascending=True)


    ##Proximity
    plt.figure()
    plt.axhline(y=0,linestyle="--",linewidth=2,color="grey")
    plt.plot(vdata['disease'], vdata['zscore'], marker='o', linestyle='-',markersize=12,linewidth=2)
    plt.xticks(rotation=45,ha='right')
    plt.ylabel('Proximity Z-Score')
    plt.show()


def plot_amspl(amspl):
    metrics = ["Shortest Path","Random Walks","Biased Random Walks", "Communicability"]
    df = amspl["Shortest Path"].T
    df = df.reset_index()
    df.columns = ["disease","Shortest Path"]

    for i in range(1,len(metrics)):
        dnew = amspl[metrics[i]].T
        dnew = dnew.reset_index()
        dnew.columns = ["disease",metrics[i]]

        df = pd.merge(df,dnew,on="disease")

    df = df.sort_values(by="Shortest Path")


    plt.figure()
    for m in metrics:
        plt.plot(df['disease'], df[m]/max(df[m]), marker='o',
                 linestyle='-',markersize=12,linewidth=2,label=m)

    plt.ylabel("Normalized AMSPL")
    plt.xticks(rotation=45,ha='right')
    plt.legend(loc="best",frameon=False)
    plt.show()


    corr = df.drop('disease', axis=1).astype(float).corr(method='spearman')
    plt.figure(figsize=(10, 8))

    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the correct aspect ratio
    sns.heatmap(corr, cmap=cmap, center=0, vmax=1,vmin=-1,
                annot=True, fmt=".2f", annot_kws={'size': 20,'color':'black'},
                square=True, linewidths=.5, cbar_kws={"shrink": 1})

    # Rotate x-axis labels and set them
    plt.xticks(np.arange(len(corr.columns)) + 0.5, corr.columns, rotation=45, ha='right')
    plt.yticks(np.arange(len(corr.columns)) + 0.5, corr.columns, rotation=0)

    # Show the plot
    plt.tight_layout()  # Adjusts plot to ensure everything fits without overlap
    plt.show()


def save(obj, file):
    with open('output/' + file,"wb") as file:
        pickle.dump(obj,file)

def load(file):
    with open('output/' + file,"rb") as file:
        obj = pickle.load(file)

    return obj

if __name__=="__main__":

    ## 1) LOAD DATA

    #Load PPI network
    with open("data/ppi_network.pkl","rb") as file:
      ppi = pickle.load(file)

    #Load drug targets
    with open("data/vitd_targets.pkl","rb") as file:
      targets = pickle.load(file)

    #Load disease genes
    with open("data/disease_genes.pkl","rb") as file:
      disease_genes = pickle.load(file)



    ## 2) EXTRACT AND EVALUATE DISEASE MODULES
    lcc_size = pd.DataFrame(columns = ["disease","size","zscore","pval"])

    for d,genes in disease_genes.items():
        data = netmedpy.lcc_significance(ppi, genes,
                                         null_model="degree_match",n_iter=10000)

        new_line = [d,data["lcc_size"],data["z_score"],data["p_val"]]
        lcc_size.loc[len(lcc_size.index)] = new_line

    save(lcc_size,"lcc_size.pkl")
    #lcc_size = load("lcc_size.pkl")

    plot_lcc_significance(lcc_size)

    #Keep only diseases with an LCC larger than 10 and statistically significant
    significant = lcc_size.query("size > 10 and zscore > 2 and pval<0.05")
    disease_names = significant.disease


    dgenes = {}
    for n in disease_names:
        lcc = netmedpy.extract_lcc(disease_genes[n], ppi)
        dgenes[n] = set(lcc.nodes)


    ## 3) CALCULATE SHORTEST PATH DISTANCE BETWEEN NODE PAIRS
    sp_distance = netmedpy.all_pair_distances(ppi,distance="shortest_path",
                                              n_processors=20,n_tasks=2000)


    ## 4) EVALUATE AMSPL BETWEEN INFLAMMATION AND FACTOR IX DEFICIENCY DISEASE
    inflammation = netmedpy.proximity(ppi, targets,
                                      dgenes["Inflammation"], sp_distance,
                                      null_model="degree_match",n_iter=10000,
                                      symmetric=False)

    factorix = netmedpy.proximity(ppi, targets,
                                      dgenes["Factor IX Deficiency"], sp_distance,
                                      null_model="degree_match",n_iter=10000,
                                      symmetric=False)

    save( (inflammation,factorix),"inf_fix.pkl" )
    #(inflamation,factorix) = load("inf_fix.pkl")

    plot_histograms(inflammation, factorix)


    ### 5) CALCULATE PROXIMITY FROM VITAMIN D TO ALL DISEASES
    vit_d = {"Vitamin D":targets}
    screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                     sp_distance,score="proximity",
                                     properties=["z_score","raw_amspl"],
                                     null_model="degree_match",
                                     n_iter=10000,n_procs=20)

    save(screen_data,"screen.pkl")
    #screen_data = load("screen.pkl")
    plot_screening(screen_data['z_score'])


    # 5) EVALUATE AMSPL UNDER DIFFERENT DISTANCES
    #Shortest Paths
    amspl = {"Shortest Path":screen_data["raw_amspl"]}

    #Random Walks
    sp_distance = netmedpy.all_pair_distances(ppi,distance="random_walk")
    screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                     sp_distance,score="proximity",
                                     properties=["raw_amspl"],
                                     null_model="degree_match",
                                     n_iter=10,n_procs=20)

    amspl["Random Walks"] = screen_data["raw_amspl"]

    #Biased Random Walks
    sp_distance = netmedpy.all_pair_distances(ppi,distance="biased_random_walk")
    screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                     sp_distance,score="proximity",
                                     properties=["raw_amspl"],
                                     null_model="degree_match",
                                     n_iter=10,n_procs=20)

    amspl["Biased Random Walks"] = screen_data["raw_amspl"]


    #Communicability
    sp_distance = netmedpy.all_pair_distances(ppi,distance="communicability")
    screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                     sp_distance,score="proximity",
                                     properties=["raw_amspl"],
                                     null_model="degree_match",
                                     n_iter=10,n_procs=20)

    amspl["Communicability"] = screen_data["raw_amspl"]


    save(amspl,"amspl.pkl")
    #amspl = load("amspl.pkl")
    plot_amspl(amspl)
