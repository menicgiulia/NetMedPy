# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 15:47:40 2024

@author: aalda
"""

import networkx as nx
import numpy as np
import pandas as pd

import NetMedPy as netmedpy
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
    plt.figure(figsize=(11,12))
    plt.bar(df['disease'], df['size'])
    plt.axhline(y=10,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC Size')
    plt.xticks(rotation=90, ha='right',va='top',fontsize=20)  # Rotate labels 45 degrees and align them right
    plt.yscale("log")
    plt.show()

    #plot z-score
    plt.figure(figsize=(11,12))
    plt.bar(df['disease'], df['zscore'])
    plt.axhline(y=2,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC Z-Score')
    plt.xticks(rotation=90, ha='right',va='top',fontsize=20)  # Rotate labels 45 degrees and align them right
    plt.show()

    #plot p-value
    plt.figure(figsize=(11,12))
    plt.bar(df['disease'], df['pval'])
    plt.axhline(y=0.05,linestyle="--",color="grey",linewidth=2)
    plt.ylabel('LCC p-value')
    plt.xticks(rotation=90, ha='right',va='top',fontsize=20)  # Rotate labels 45 degrees and align them right
    plt.ylim(10**-5,1)
    plt.yscale("log")
    plt.show()


def plot_histograms(inflammation, huntington):

    plt.figure()
    sns.kdeplot(inflammation['dist'], color='blue', fill=True,
                alpha=0.5, label='Inflammation')
    sns.kdeplot(huntington['dist'], color='red', fill=True,
                alpha=0.5, label='Huntington')


    plt.axvline(x=inflammation['raw_amspl'],linewidth=10,
                color='blue',alpha=0.5)
    plt.text(x=inflammation['raw_amspl'], y=7,
             s=f"  Z = {inflammation['z_score']:.2f}",
             color='blue', verticalalignment='top',
             horizontalalignment='left',fontsize=30)


    plt.axvline(x=huntington['raw_amspl'],linewidth=10,color='red',alpha=0.5)

    plt.text(x=huntington['raw_amspl'], y=7,
             s=f"  Z = {huntington['z_score']:.2f}",
             color='red', verticalalignment='top',
             horizontalalignment='left',fontsize=30)

    plt.xlabel('AMSPL from Vitamin D')
    plt.ylabel('Density')
    plt.ylim(0,9)
    plt.legend(loc='best',frameon=False)
    plt.show()


def plot_screening(df):
    ##Proximity
    plt.figure(figsize=(11,12))
    plt.axhline(y=0,linestyle="--",linewidth=2,color="grey")
    plt.plot(df['Disease'], df['z_score'], marker='o', linestyle='-',markersize=12,linewidth=2)
    plt.xticks(rotation=90,ha='right',fontsize=22)
    plt.ylabel('Proximity Z-Score')
    plt.show()


def plot_amspl(df):

    metrics = df.columns
    metrics = df.columns[2:len(metrics)]

    plt.figure(figsize=(11,12))

    for m in metrics:
        plt.plot(df['Disease'], df[m]/max(df[m]), marker='o',
                 linestyle='-',markersize=12,linewidth=2,label=m)

    plt.ylabel("Normalized AMSPL")
    plt.xticks(rotation=90,ha='right',fontsize=22)
    plt.legend(loc="best",frameon=False)
    plt.show()


    cordf = df[list(metrics)]
    corr = cordf.astype(float).corr(method='spearman')
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
    with open('examples/VitaminD/output/' + file,"wb") as file:
        pickle.dump(obj,file)

def load(file):
    with open('examples/VitaminD/output/' + file,"rb") as file:
        obj = pickle.load(file)

    return obj


if __name__=="__main__":

    lcc_size = load("lcc_size.pkl")
    plot_lcc_significance(lcc_size)

    screen = load('screen.pkl')
    amspl = load('amspl.pkl')

    df = screen['z_score'].T
    df = df.reset_index()
    df.columns= ['Disease','z_score']


    for k,v in amspl.items():
        nd = v.T.reset_index()
        nd.columns = ['Disease',k]

        df = pd.merge(df,nd,on='Disease')

    df = df.sort_values(by='z_score',ascending=True)

    plot_screening(df)
    plot_amspl(df)
