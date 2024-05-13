# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:50:42 2023

@author: aalda
"""

import proximity.distances as old_distances
import netmedpy.NetMedPy as new_distances
import other_distances.ExperimentalDistances as expDistances
import networkx as nx
import numpy as np
import tools.LoadContext as LoadContext
import tools.Cronometer as Cronometer
import pandas as pd
import matplotlib.pyplot as plt
import random
from proximity.network import Network as agNetwork


def exe_time(CONTEXT,s,reps,method,distMatrix = None):
    ppi = agNetwork(CONTEXT.PPI) if method=='oo_dijkstra_halt' else CONTEXT.PPI
    nodes = set(CONTEXT.PPI.nodes)

    times = []
    c = Cronometer.Cronometer()

    for i in range(reps):
        A = random.sample(nodes, s)
        B = random.sample(nodes, s)

        if method == 'oo_dijkstra_halt':
            c.tick()
            prox = old_distances.get_min_shortest_paths(ppi, A, B)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        elif method == 'oa_dijkstra_no_halt': # method == 'new'
            c.tick()
            prox = expDistances.get_avg_min_shortest_path_oa_nohalt(ppi, A, B)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")

        elif method == 'oa_dijkstra_halt':
            c.tick()
            prox = expDistances.get_avg_min_shortest_path_oa_halt(ppi,A,B)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                   ",end="")
        elif method == 'oo_dijkstra_halt_2':
            c.tick()
            prox = new_distances._get_amspl_all_operations(ppi, A, B)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        elif method == 'distance_matrix':
            c.tick()
            prox = new_distances.get_amspl(ppi, A, B,distMatrix)
            c.tock()
            print(f"\rSize {s} rep {i} time {c.elapsed_seconds}                    ",end="")
        else:
            raise ValueError(f'Method {method} not recognized')

        times.append(c.elapsed_seconds)
    return np.array(times).mean()




def create_serie(CONTEXT,size_ini,size_fin,size_step, reps = 10,method='old',dmatrix=None):
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

    outfile = f"C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/proximity/singleCallSeries/{method}_{reps}.csv"

    serie.to_csv(outfile,index=False)



def plot_series(reps):
    inFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/proximity/singleCallSeries/"

    oa_d_h = pd.read_csv(inFile + f"oa_dijkstra_halt_{reps}.csv")
    oa_d_nh = pd.read_csv(inFile + f"oa_dijkstra_no_halt_{reps}.csv")
    oo_d_h = pd.read_csv(inFile + f"oo_dijkstra_halt_{reps}.csv")
    oo_d_h_2 = pd.read_csv(inFile + f"oo_dijkstra_halt_2_{reps}.csv")
    dmat = pd.read_csv(inFile + f"distance_matrix_{reps}.csv")

    plt.figure()
    plt.plot(oo_d_h['size'],oo_d_h['time'],label='One to One - halt',color = 'black',
             linestyle='-',linewidth = 2)
    # plt.plot(oo_d_h_2['size'],oo_d_h_2['time'],label='One to One - halt v.2',color = 'black',
    #          linestyle='--',linewidth = 1)
    # plt.plot(oa_d_h['size'],oa_d_h['time'],label='One to all - halt')
    plt.plot(oa_d_nh['size'],oa_d_nh['time'],label='One to all - No halt',color = 'blue')
    plt.plot(dmat['size'],dmat['time'],label='Distance Matrix',color="red",linewidth=2)
    plt.xlabel("A,B Size")
    plt.ylabel(f"Execution Time (s) ({reps} reps)")
    plt.legend(loc="best",frameon=False)
    plt.xlim(left=0)
    plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    #CONTEXT = LoadContext.Context()
    size_ini = 100
    size_fin = 2001
    size_step = 100
    reps = 5

    # print('oo_dijkstra_halt') #Original Rodrigo implementation
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                              'oo_dijkstra_halt', reps)

    # print('oa_dijkstra_no_halt')
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                              'oa_dijkstra_no_halt', reps)

    # print('oa_dijkstra_halt')
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                              'oa_dijkstra_halt', reps)


    # print('oo_dijkstra_halt_2') #Own implementation
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'oo_dijkstra_halt_2', reps)

    # print('dist_matrix')
    # dmat = new_distances.load_distances(CONTEXT.DISTANCE_FILE)
    # create_execution_time_series(CONTEXT, size_ini, size_fin, size_step,
    #                               'distance_matrix', reps,dmatrix = dmat)




    plot_series(reps)
