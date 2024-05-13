# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 18:29:34 2023

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


if __name__=='__main__':
    CONTEXT = LoadContext.Context()

    reps = 10

    c = Cronometer.Cronometer()

    times_calc = []

    for i in range(reps):

        c.tick()
        distances = new_distances.all_pair_distances(CONTEXT.PPI)
        c.tock()
        secs = c.elapsed_seconds
        times_calc.append(secs)
        print(f"Calc {c.format_seconds()}")

    avg_calc = np.array(times_calc).mean()

    print(f"Avg Calc time {avg_calc}")


    c.tick()
    new_distances.save_distances(distances, CONTEXT.DISTANCE_FILE)
    c.tock()
    print(f"Write {c.format_seconds()}")


    times_read = []

    for i in range(reps):

        c.tick()
        distances = new_distances.load_distances(CONTEXT.DISTANCE_FILE)
        c.tock()
        secs = c.elapsed_seconds
        times_read.append(secs)
        print(f"read {c.format_seconds()}")
    avg_read = np.array(times_read).mean()

    print(f"Avg read time {avg_read}")

    df = pd.DataFrame({"Calc":[avg_calc],
                       "Load":[avg_read]
        })

    file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/Postdoc/proximitySeparationOptimization/proximity/data/all_pair_distances/avg_time.csv"

    df.to_csv(file,index=False)
