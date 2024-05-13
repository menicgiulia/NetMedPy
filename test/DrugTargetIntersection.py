# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 00:47:53 2024

@author: aalda
"""

import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt

dbFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/Drug_Bank/DB_Drug_Targets_2023.csv"

db = pd.read_csv(dbFile)
db = db.query("organism=='Humans'")
db = db[db['Status'].str.contains('approved',case=False)]

names = list(db.DB_id.unique())

tries = 30000
hits = 0
ints = []

for i in range(tries):
    a = random.choice(names)
    b = random.choice(names)

    ta = set(db.query("DB_id==@a")['Gene_Target'])
    tb = set(db.query("DB_id==@b")['Gene_Target'])

    il = len(ta & tb)
    if il != 0:
        hits = hits+1
        ints.append(il)
        
prob = hits/tries

print(f"Intersection probability {prob:.2f}")

plt.figure(figsize=(9,6))
plt.hist(ints,bins=25)
plt.tick_params(axis='both', which='major', labelsize=14) 
plt.ylabel("Frequency",fontsize=14)
plt.xlabel("Shared Genes",fontsize=14)
plt.yscale("log")
plt.show()

npa = np.array(ints)

np.median(npa)


genes = list(db.Gene_Target.unique())

distGenes = []
for g in genes:
    l = db.query("Gene_Target==@g")
    d = len(l.DB_id.unique())
    distGenes.append(d)
    
df = pd.DataFrame({'Gene':genes,'A_drugs':distGenes})
    
plt.figure(figsize=(9,6))
plt.hist(distGenes,bins=50)
plt.tick_params(axis='both', which='major', labelsize=14) 
plt.ylabel("Frequency",fontsize=14)
plt.xlabel("Associated Drugs",fontsize=14)
plt.yscale("log")
plt.show()