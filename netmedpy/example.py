# -*- coding: utf-8 -*-
"""
Example of Functionality and Application of NetMedPy

This example demonstrates how to utilize the NetMedPy package. Ensure the following 
dependencies are installed before proceeding:

Required packages:
- networkx
- pandas
- numpy
- pickle
- multiprocessing
- random
- scipy
- ray
- seaborn

Setup and Execution Instructions:

1) Installation:

   NetMedPy is compatible with Python > 3.6 and < 3.12

   If NetMedPy was installed via pip, skip to step 2. Otherwise, manually verify that you have 
   installed the packages mentioned above. Add the NetMedPy directory 
   to your PYTHONPATH environment variable to ensure Python can locate the package:

   Linux/Mac:

    
    1) If you installed NetMedPy using pip, ignore this step. Otherwise
        include the netmedpy directory in the PYTHONPATH variable:
            
        On Linux/Mac:
            
            export PYTHONPATH="/user_path_to/NetMedPy/netmedpy:$PYTHONPATH"
        
        On Windows:
        
            set PYTHONPATH=C:\user_path_to\NetMedPy\netmedpy;%PYTHONPATH%

    2) Navigate to the NetMedPy directory where this example script is located.

    3) Run the script using Python > 3.6 and < 3.12:
            python example.py
            
    4) Enjoy

For any issues or further instructions, refer to the documentation or contact the authors.
    
@author: Andres Aldana Gonzalez (a.aldana@northeastern.edu)
"""
import networkx as nx
import NetMedPy as netmedpy
import random
import seaborn as sns
import matplotlib.pyplot as plt
from time import sleep


def plot_proximity_distribution(p):
    plt.figure(figsize=(9, 6))
    sns.kdeplot(p['dist'], fill=True, label='KDE of ASPL')
    plt.axvline(x=p['raw_amspl'], color="red", linewidth=2, label="AMSPL(S,T)")
    plt.text(p['raw_amspl'], 0, f"z-score:{p['z_score']:.2f} p:{p['p_value_single_tail']:.2f}",
             color="black", ha='right',
             verticalalignment="bottom", transform=plt.gca().get_xaxis_transform())
    plt.xlabel("ASPL")
    plt.suptitle("Proximity")
    plt.legend(loc="best", frameon=False)
    plt.show()


def plot_separation_distribution(s):
    plt.figure(figsize=(9, 6))
    sns.kdeplot(s['dist'], fill=True, label='KDE of Separation')
    plt.axvline(x=s['raw_separation'], color="red", linewidth=2, label="Separation(S,T)")
    plt.text(s['raw_separation'], 0, f"z-score: {s['z_score']:.2f} p:{s['p_value_double_tail']:.2f}",
             color="black", ha='right',
             verticalalignment="bottom", transform=plt.gca().get_xaxis_transform())
    plt.xlabel("Separation")
    plt.suptitle("Separation")
    plt.legend(loc="best", frameon=False)
    plt.show()


def plot_lcc_distribution(lcc_data):
    plt.figure(figsize=(9, 6))
    sns.kdeplot(lcc_data['dist'], fill=True, label='KDE of LCC size')
    plt.axvline(x=lcc_data['lcc_size'], color="red", linewidth=2, label="G1 LCC size")
    plt.text(lcc_data['lcc_size'], 0, f"z-score: {lcc_data['z_score']:.2f}", color="black", ha='right',
             verticalalignment="bottom", transform=plt.gca().get_xaxis_transform())
    plt.xlabel("LCC Size")
    plt.suptitle("LCC Significance")
    plt.legend(loc="best", frameon=False)
    plt.show()


if __name__ == "__main__":

    # --- NETWORK CREATION ---
    m = 3  # Number of edges to attach from a new node to existing nodes
    G1 = nx.barabasi_albert_graph(300, m)
    G2 = nx.barabasi_albert_graph(700, m)
    mapping = {node: node + 300 for node in G2.nodes()}
    G2 = nx.relabel_nodes(G2, mapping)
    for _ in range(10):
        node_from_G1 = random.choice(list(G1.nodes()))
        node_from_G2 = random.choice(list(G2.nodes()))
        G1.add_edge(node_from_G1, node_from_G2)
    G = nx.compose(G1, G2)

    # --- SAMPLING NODES ---
    S = random.sample(list(G1.nodes()), 20)
    T = random.sample(list(G2.nodes()), 20)

    # --- CALCULATING DISTANCES ---
    D = netmedpy.all_pair_distances(G, distance="shortest_path", n_processors=10, n_tasks=100)
    mat = D.matrix

    
    netmedpy.save_distances(D, "distances.pkl") # If you want to save the distance matrix
    D = netmedpy.load_distances("distances.pkl") # If you want to load a pre-computed distance matrix

    # --- PROXIMITY ANALYSIS ---
    p = netmedpy.proximity(G, S, T, D, null_model='degree_match', n_iter=1000)
    plot_proximity_distribution(p)

    # --- SEPARATION ANALYSIS ---
    rawSep = netmedpy.separation(G, S, T, D)
    print(f"Separation: {rawSep:2f}")
    s = netmedpy.separation_z_score(G, S, T, D, null_model='degree_match', n_iter=1000)
    plot_separation_distribution(s)

    # --- LCC ANALYSIS ---
    L = list(G1.nodes())
    lcc = netmedpy.extract_lcc(L, G)
    lcc_data = netmedpy.lcc_significance(G, L, null_model='degree_match', n_iter=10000)
    print(f"LCC-size={lcc_data['lcc_size']} z-score={lcc_data['z_score']:0.2f} p-value={lcc_data['p_val']:0.2f}")
    plot_lcc_distribution(lcc_data)



    # --- SELECTING DIFFERENT METRICS ---
    D = netmedpy.all_pair_distances(G, distance="random_walk",reset=0.1)
    print("Distance: Random Walk")
    print(D.matrix[:5,:5])

    D = netmedpy.all_pair_distances(G, distance="communicability")
    print("Distance: communicability")
    print(D.matrix[:5,:5])


    #User defined distance.
    #Users can define their own distances by setting distance='custom' and custom_distance= function that calculates the distance between a single node and every other node.
    #The result should be a dictionary with the distance between the source node and all other nodes in the network.
    def d(source, graph, v):
        d = {}

        for b in graph.nodes():
            d[b] = v+1

        return d

    D = netmedpy.all_pair_distances(G, distance='custom',custom_distance=d,v=3)

    print("Distance: user defined")
    print(D.matrix[:5,:5])
    
    sleep(30) # Gives ray time to shutdown properly
    
    # --- SCREEN SOURCE AND TARGET NODES ---
    
    #Define source nodes
    sources = {}
    sources["S1"] = random.sample(list(G1.nodes()), 20)
    sources["S2"] = random.sample(list(G2.nodes()), 20)
    
    #Define target nodes
    targets = {}
    targets["T1"] = random.sample(list(G1.nodes()), 20)
    targets["T2"] = random.sample(list(G2.nodes()), 20)
    
    
    #Calculate distance matrix
    D = netmedpy.all_pair_distances(G, distance="shortest_path",n_processors=10,n_tasks=100)
    
    sleep(30)

    #Calculate proximity z-score between sources and targets
    screening = netmedpy.screening(sources, targets, G, D,score="proximity",properties=["z_score"],
                       null_model="degree_match",n_iter=2000,symmetric=True,
                       n_procs=2)
    print("Screening with proximity")
    print("Z-Score")
    print(screening["z_score"])

    sleep(30)

    #Separation z_score using uniform sampling and retrieving different scores
    screening = netmedpy.screening(sources, targets, G, D,score="separation_z_score",
                       properties=["z_score","raw_separation","p_value_single_tail"],
                       null_model="uniform",n_iter=2000,symmetric=False,
                       n_procs=2)
    
    print("Screening with separation")
    print("Z-Score")
    print(screening["z_score"])
    print("Separation")
    print(screening["raw_separation"])
    print("Single tail p-value")
    print(screening["p_value_single_tail"])
