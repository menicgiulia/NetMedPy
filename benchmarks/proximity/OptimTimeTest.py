# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 16:53:27 2023

@author: aalda
"""

import proximity.distances as old_distances
import netmedpy.NetMedPy as new_distances
import networkx as nx
import numpy as np
import tools.LoadContext as LoadContext
import tools.Cronometer as Cronometer
import pandas as pd
import matplotlib.pyplot as plt
import random
from proximity.network import Network as agNetwork
import tools.General as genTools
from multiprocessing import Pool, cpu_count
from tqdm import tqdm



def exe_time(CONTEXT,s,reps,method,distMatrix = None):
    ppi = CONTEXT.PPI
    nodes = set(CONTEXT.PPI.nodes)

    times = []
    c = Cronometer.Cronometer()

    for i in range(reps):
        A = random.sample(nodes, s)
        B = random.sample(nodes, s)

        if method == 'old_proximity':
            c.tick()
            prox = new_distances.proximity_allCalc(ppi, A, B,'exact',1000,0)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        elif method == 'distance_matrix':
            c.tick()
            prox = new_distances._proximity_dmatrix(ppi, A, B, distMatrix,'exact',1000,0)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")

        elif method == 'parallel':
            c.tick()
            prox = new_distances.proximity_dm_parallel(ppi, A, B, distMatrix, -1,n_iter=1000)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                   ",end="")

        else:
            raise ValueError(f"Method {method} not recognized")

        times.append(c.elapsed_seconds)
    return np.array(times).mean()




def create_serie(CONTEXT,size_ini,size_fin,size_step, reps = 10,method=None,dmatrix=None):
    sizes = np.arange(size_ini,size_fin,size_step)
    times = np.zeros(len(sizes))


    for i in range(len(sizes)):
        s = sizes[i]
        times[i] = exe_time(CONTEXT,s,reps,method,distMatrix=dmatrix)
    print("")

    res = pd.DataFrame({'size':sizes,'time':times})
    return res


def create_execution_time_series(CONTEXT,size_ini,size_fin,size_step,method,reps,dmatrix=None):
    serie = create_serie(CONTEXT,size_ini,size_fin,size_step, reps = reps,method=method,dmatrix=dmatrix)

    outfile = f"C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/proximity/OptimTime/{method}_{reps}.csv"

    serie.to_csv(outfile,index=False)



def plot_series(reps):
    inFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/proximity/OptimTime/"

    old = pd.read_csv(inFile + "old_proximity_1.csv")
    dm = pd.read_csv(inFile + f"distance_matrix_{reps}.csv")


    plt.figure()
    plt.plot(old['size'],old['time'],label='All calculations - 1 rep',color = 'black',
              linestyle='-',linewidth = 2)
    plt.plot(dm['size'],dm['time']+35,label=f"Distance Matrix - {reps} reps",color = 'red', linewidth=2)
    plt.xlabel("A,B Size")
    plt.ylabel("Execution Time (s)")
    plt.legend(loc="best",frameon=False)
    # plt.yscale("log")
    plt.xlim(left=0)
    plt.show()


if __name__ == '__main__':
    CONTEXT = LoadContext.Context()
    reps = 5

    size_ini = 50
    size_fin = 300
    size_step = 50

    # c = Cronometer.Cronometer()

    # print('old proximity') #Original Rodrigo implementation
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'old_proximity', 1,None)

    # c.tick()
    # distances = new_distances.load_distances(CONTEXT.DISTANCE_FILE)
    # c.tock()
    # print(f"Elapsed {c.format_seconds()}")

    # print('distance matrix')
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'distance_matrix', reps,distances)

    # print('parallel')
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                                'parallel', reps,distances)


    plot_series(reps)
