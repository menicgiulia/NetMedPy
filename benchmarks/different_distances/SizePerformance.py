# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 04:42:01 2023

@author: aalda
"""
import networkx as nx
import netmedpy.NetMedPy as network_metrics
import netmedpy.DistanceMatrix as DistanceMatrix
import numpy as np
import tools.Cronometer as Cronometer
import pandas as pd


import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-colorblind')
plt.rcParams["figure.figsize"] = [9, 6]
plt.rcParams["figure.autolayout"] = True

plt.rcParams['font.size'] = 15
plt.rcParams.update({'font.size':15})
# Set the axes labels font size
plt.rc('axes', labelsize=15)
# Set the font size for x tick labels
plt.rc('xtick', labelsize=15)
# Set the font size for y tick labels
plt.rc('ytick', labelsize=15)



def calc_row(size, iterations):
    m = 2
    G_ba = nx.barabasi_albert_graph(size, m, seed=42)

    row = [size]
    methods = ["shortest_path","random_walk","biased_random_walk","communicability"]

    for m in methods:
        times = []
        c = Cronometer.Cronometer()

        for j in range(iterations):

            print(f"Size={size} Iteration={j+1} Method={m}")
            c.tick()
            D = network_metrics.all_pair_distances(G_ba, distance=m,reset = 0.2)
            c.tock()

            times.append(c.elapsed_milliseconds())
            del D

        t = np.mean(times)
        row.append(t)

    return row


def plot_data(file):
    df = pd.read_csv(file)

    plt.figure()
    plt.plot(df.Size,df.shortest_path/1000,linewidth=2,label="Shortest Path",
             color="black")
    plt.scatter(df.Size,df.shortest_path/1000,s=80,marker="o",color="black")

    plt.plot(df.Size,df.random_walk/1000,linewidth=2,label="Random  Walk",
             color="blue")
    plt.scatter(df.Size,df.random_walk/1000,s=80,marker="o",color="blue")

    plt.plot(df.Size,df.biased_random_walk/1000,linewidth=2,
             label="Biased Random Walk",color="purple")
    plt.scatter(df.Size,df.biased_random_walk/1000,s=80,marker="o",color="purple")

    plt.plot(df.Size,df.communicability/1000,linewidth=2,label="Communicability",
             color="darkgreen")
    plt.scatter(df.Size,df.communicability/1000,s=80,marker="o",color="darkgreen")

    plt.xlabel("Size")
    plt.ylabel("Time(s)")
    plt.legend(loc="best")
    plt.show()


if __name__=="__main__":

    output_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/all_pair_distances/size_time_2.csv"

    ini_size = 1000
    increment = 1000
    final_size = 18000
    iterations = 3

    df = pd.DataFrame(columns=["Size","shortest_path","random_walk",
                                "biased_random_walk","communicability"])

    size = ini_size
    while(size <= final_size):
        row = calc_row(size,iterations)

        df.loc[len(df.index)] = row
        size += increment

    df.to_csv(output_file,index=False)
    print("")

    plot_data(output_file)
