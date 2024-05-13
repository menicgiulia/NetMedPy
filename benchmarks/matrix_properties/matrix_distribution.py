# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:37:28 2024

@author: aalda
"""

import netmedpy.NetMedPy as netmedpy
import netmedpy.DistanceMatrix as dmatrix
import matplotlib.pyplot as plt

import random

import numpy as np

mat = "D:/data/distanceMatrix/channing/ppi932023/communicability.pkl"

D = netmedpy.load_distances(mat)

m = D.matrix


array = m.flatten()

sample_size = 10000

random_sample = np.random.choice(array, size=sample_size, replace=False)

# random_sample = np.log(random_sample[random_sample < 20])

# random_sample = np.log(random_sample)

plt.figure(figsize=(9,6))
plt.hist(random_sample, bins=50,density=True)
plt.ylabel("Density")
plt.xlabel("log(Distance)")
plt.show()

random_sample.mean()
np.median(random_sample)
random_sample.std()

random_sample.max()
