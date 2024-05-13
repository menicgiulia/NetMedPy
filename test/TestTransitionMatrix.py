import numpy as np
import networkx as nx
from numpy.linalg import inv
import tools.Cronometer as Cronometer
import pickle

import matplotlib.pyplot as plt

import netmedpy.DistanceMatrix as dMatrix
import netmedpy.NetMedPy as netmedpy

def RWR_Matrix(G, c):
    # Map from node names to indices and vice versa
    node_to_index = {node: i for i, node in enumerate(G.nodes())}
    index_to_node = {i: node for node, i in node_to_index.items()}

    A = nx.to_numpy_array(G)

    # Normalize A to get transition probabilities
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Avoid division by zero for isolated nodes
    T = A / col_sums

    # Calculate the modified transition matrix for RWR
    n = G.number_of_nodes()
    I = np.eye(n)
    M = c * np.linalg.inv(I - (1 - c) * T)
    M = M.T

    return (M,node_to_index,index_to_node)


def random_walk_with_restart(start_node, RWR_data):

    M,node_to_index,index_to_node = RWR_data
    # Initial distribution vector
    n = len(node_to_index.keys())
    r0 = np.zeros(n)
    r0[node_to_index[start_node]] = 1  # Set starting node

    # Calculate the steady-state distribution r
    r = M.dot(r0)

    # Normalize r to ensure it sums to 1
    r /= np.sum(r)

    # Map the result back to node names
    r_named = {index_to_node[i]: val for i, val in enumerate(r)}

    return r_named


# Example usage
if __name__ == "__main__":
    # Create a weighted graph with named nodes
    # G = nx.Graph()
    # G.add_edge('A', 'B', weight=1)
    # G.add_edge('B', 'C', weight=1)
    # G.add_edge('C', 'D', weight=1)
    # G.add_edge('D', 'A', weight=1)
    # G.add_edge('A', 'C', weight=2)  # Example of a weighted edge
    # G.add_edge('C','E',weight=0.1)
    # G.add_edge('B','F',weight=10)

    G = nx.barabasi_albert_graph(20, 2)


    # Name of the starting node
    start_node = 0

    nx.draw(G,with_labels=True)

    # Restart probability
    c = 0.2


    # Calculate the RWR
    d = RWR_Matrix(G, c)

    r_named = random_walk_with_restart(start_node, d)


    (M,node_to_index,index_to_node) = d

    print("Modified Transition Matrix M:/n", M)
    print("Steady State Vector r (by node names):/n", r_named)

    print(f"The sum {sum(r_named.values())}")



    f = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"
    with open(f,'rb') as file:
        ppi = pickle.load(file)


    ppi = netmedpy.extract_lcc(ppi.nodes(), ppi)

    c = Cronometer.Cronometer()

    c.tick()
    (M,node_to_index,index_to_node) = RWR_Matrix(ppi, 0.25)
    c.tock()

    print(f"Time: {c.format_seconds()}")


    # normalized = np.log(1 + 1/(1+M))

    # normalized = np.log(1 + 1/(0.0000000001+M**2))
    # normalized = normalized/np.max(normalized)

    normalized = np.log(M)

    ma = np.max(normalized)
    mi = np.min(normalized)

    normalized = 1 - ((normalized - mi)/(ma - mi))

    array = normalized.flatten()

    sample_size = 2000

    random_sample = np.random.choice(array, size=sample_size, replace=False)

    plt.figure(figsize=(9,6))
    plt.hist(random_sample, bins=50,density=True)
    plt.ylabel("Density")
    plt.xlabel("Distance")
    plt.xlim(0,1)
    plt.show()


    dsorted = np.sort(random_sample)
    cumProb = np.arange(1,len(random_sample)+1)/len(random_sample)

    plt.figure(figsize=(9,6))
    plt.step(dsorted, cumProb,where='mid', lw=2)
    plt.ylabel("CDF")
    plt.xlabel("Distance")
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.show()



    #####################For Biased Random Walks################

    G = ppi
    degs = []
    for i in range(len(G)):
        degs.append(G.degree(index_to_node[i]))

    degs = np.array(degs)

    B = M / degs
    B = B.T
    B = B / B.sum(axis=0)
    B = B.T


    #normalized = np.log(1 + 1/(0.0000000001+B**2))
    #normalized = normalized/np.max(normalized)

    normalized = np.log(B)

    ma = np.max(normalized)
    mi = np.min(normalized)

    normalized = 1 - ((normalized - mi)/(ma - mi))

    array = normalized.flatten()

    sample_size = 2000

    random_sample = np.random.choice(array, size=sample_size, replace=False)

    plt.figure(figsize=(9,6))
    plt.hist(random_sample, bins=50,density=True)
    plt.ylabel("Density")
    plt.xlabel("Distance")
    plt.xlim(0,1)
    plt.show()


    dsorted = np.sort(random_sample)
    cumProb = np.arange(1,len(random_sample)+1)/len(random_sample)

    plt.figure(figsize=(9,6))
    plt.step(dsorted, cumProb,where='mid', lw=2)
    plt.ylabel("CDF")
    plt.xlabel("Distance")
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.show()



    ###################For Comunicability


    D = dMatrix.DistanceMatrix()
    D._from_name_list(list(ppi.nodes()))
    comm = nx.communicability_exp(ppi)

    for a,d in comm.items():
        for b,c in d.items():
            D.put(a,b,c)
    del comm

    # D.matrix = 1/D.matrix

    # median = np.median(D.matrix)
    # D.matrix = D.matrix/median


    normalized = D.matrix

    normalized = np.log(normalized)

    ma = np.max(normalized)
    mi = np.min(normalized)

    normalized = 1 - ((normalized - mi)/(ma - mi))


    array = normalized.flatten()

    sample_size = 2000

    random_sample = np.random.choice(array, size=sample_size, replace=False)

    plt.figure(figsize=(9,6))
    plt.hist(random_sample, bins=50,density=True)
    plt.ylabel("Density")
    plt.xlabel("Distance")
    plt.xlim(0,1)
    plt.show()


    dsorted = np.sort(random_sample)
    cumProb = np.arange(1,len(random_sample)+1)/len(random_sample)

    plt.figure(figsize=(9,6))
    plt.step(dsorted, cumProb,where='mid', lw=2)
    plt.ylabel("CDF")
    plt.xlabel("Distance")
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.show()
