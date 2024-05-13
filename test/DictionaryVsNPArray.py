# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 00:56:45 2023

@author: aalda
"""
import numpy as np
import random
import tools.Cronometer as cronometer
import matplotlib.pyplot as plt
import pandas as pd

def create_dictionary(size):
    D = {}
    for i in range(size):
        ki = f"a{i}"

        d = {}
        for j in range(size):
            kj = f"a{j}"
            d[kj] = random.uniform(0,100)
        D[ki] = d
    return D


class DistanceMatrix:
    def __init__(self, distance_dict):
        # Extract unique nodes from the dictionary keys
        self.nodes = list(distance_dict.keys())

        # Create a mapping from node name to integer index
        self.node_to_idx = {node: idx for idx, node in enumerate(self.nodes)}

        # Initialize the 2D array with inf (since not every node will be connected to every other node)
        self.array = np.full((len(self.nodes), len(self.nodes)), np.inf)

        # Populate the array using the distance dictionary
        for node_i, neighbors in distance_dict.items():
            for node_j, distance in neighbors.items():
                self.array[self.node_to_idx[node_i], self.node_to_idx[node_j]] = distance

    def get_distance(self, A, B):
        # Retrieve the indices of nodes A and B
        idx_A = self.node_to_idx.get(A)
        idx_B = self.node_to_idx.get(B)

        # Check if nodes A and B are valid
        if idx_A is None or idx_B is None:
            raise ValueError(f"One or both nodes {A} and {B} are not present in the matrix.")

        # Return the distance from the 2D array
        return self.array[idx_A, idx_B]


def calcTimes(size,chunk):
    df = pd.DataFrame(columns=["Operations","Time_D","Time_A"])

    print("Creating dictionary")
    D = create_dictionary(size)

    print("Creating Distance Matrix")
    mat = DistanceMatrix(D)

    print("Done")

    cron = cronometer.Cronometer()

    for i in range(1,201):

        operations = i*chunk

        cron.tick()

        for _ in range(operations):
            i1 = random.randint(0, size-1)
            i2 = random.randint(0,size-1)

            k = D[f"a{i1}"][f"a{i2}"]
        cron.tock()



        td = cron.elapsed_seconds

        cron.tick()

        for _ in range(operations):

            i1 = random.randint(0, size-1)
            i2 = random.randint(0,size-1)

            k = mat.get_distance(f"a{i1}",f"a{i2}")
        cron.tock()

        ta = cron.elapsed_seconds

        row = [operations,td,ta]

        df.loc[len(df.index)] = row
    return df

def plot_result(df):
    plt.figure()
    plt.plot(df.Operations,df.Time_D,label="Time D",color="blue")
    plt.plot(df.Operations,df.Time_A,label="Time A",color="red")
    plt.legend(loc="best",frameon=False)
    plt.show()


def main():
    df = calcTimes(8000, 1000)

    plot_result(df)

if __name__=='__main__':
    main()


s = 20000
d = np.full((s,s),1)

for i in range(s):
    for j in range(s):
        d[i][j] = random.uniform(0,100)
