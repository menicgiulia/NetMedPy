# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:50:42 2023

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


def exe_time(CONTEXT,s,reps,method,distMatrix = None):
    ppi = agNetwork(CONTEXT.PPI) if method=='traditional' else CONTEXT.PPI
    nodes = set(CONTEXT.PPI.nodes)

    times = []
    c = Cronometer.Cronometer()

    for i in range(reps):
        A = random.sample(nodes, s)
        B = random.sample(nodes, s)

        if method == 'traditional':
            c.tick()
            prox = old_distances.separation(ppi, A, B)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        elif method == 'optimized': # method == 'new'
            c.tick()
            prox = new_distances.separation(ppi,A, B, distMatrix)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        else:
            raise ValueError(f'Method {method} not recognized')

        times.append(c.elapsed_seconds)
    return np.array(times).mean()




def create_serie(CONTEXT,size_ini,size_fin,size_step, reps = 10,method='traditional',dmatrix=None):
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

    outfile = f"C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/separation/singleCall/{method}_{reps}.csv"

    serie.to_csv(outfile,index=False)



def plot_series(reps):
    inFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/separation/singleCall/"

    trad = pd.read_csv(inFile + f"traditional_{reps}.csv")
    optim = pd.read_csv(inFile + f"optimized_{reps}.csv")

    plt.figure()
    plt.plot(trad['size'],trad['time'],label='Original',color = 'black',
             linestyle='-',linewidth = 2)
    plt.plot(optim['size'],optim['time'],label='Optimized',color = 'red',
             linestyle='-',linewidth = 2)

    plt.xlabel("A,B Size")
    plt.ylabel(f"Execution Time (s) ({reps} reps)")
    plt.legend(loc="best",frameon=False)
    plt.xlim(left=0)
    # plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    CONTEXT = LoadContext.Context()
    size_ini = 100
    size_fin = 1500
    size_step = 100
    reps = 5

    # print('Traditional') #Original Rodrigo implementation
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'traditional', reps)

    # print('dist_matrix')
    # dmat = new_distances.load_distances(CONTEXT.DISTANCE_FILE)
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'optimized', reps,dmatrix = dmat)


    plot_series(reps)
