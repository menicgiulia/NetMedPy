# -*- coding: utf-8 -*-
# Uncomment the following lines to set the working directory to the VitaminD directory.
# If commented, navigate to the directory manually.
# import os
# os.chdir("/user_path_to/NetMedPy/examples/VitaminD/")

### Import necessary packages
import os
import pickle
import random
import time as time_functions
import numpy as np
import pandas as pd
import ray
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as opt


from Cronometer import Cronometer as cronometer
import guney.network as p_network
import netmedpy

### Plot style
plt.style.use('seaborn-v0_8-colorblind')
plt.rcParams["figure.figsize"] = [9, 6]
plt.rcParams["figure.autolayout"] = True
font = 20

plt.rcParams['font.size'] = font
plt.rcParams.update({'font.size': font})
plt.rc('axes', labelsize=font)
plt.rc('xtick', labelsize=font)
plt.rc('ytick', labelsize=font)

# Specify the number of cores used for calculations
n_procs = 20

### Auxiliary functions for saving and loading data
def save(obj, file):
    with open(file, "wb") as file:
        pickle.dump(obj, file)


def load(file):
    with open(file, "rb") as file:
        return pickle.load(file)


def savefig(fname):
    if not os.path.exists("output/plots/"):
        os.makedirs("output/plots/")
    plt.savefig(f'output/plots/{fname}.png', format='png')


### Load PPI network, Vitamin D targets, and Disease genes:
# Load PPI network
ppi = load("data/input/ppi/ppi_network.pkl")

# Load drug targets
targets = load("data/input/drug_targets/vitd_targets.pkl")
targets = targets & set(ppi.nodes)

# Load disease genes
disease_genes = load("data/input/disease_genes/disease_genes.pkl")
for d, g in disease_genes.items():
    disease_genes[d] = set(g) & set(ppi.nodes)


### Evaluate the size and statistical significance of the largest connected component of each disease
# Function to plot the LCC statistical significance
def plot_lcc_significance(lcc_data):
    df = lcc_data.sort_values(by="size", ascending=False).fillna(0)

    # Plot size
    plt.figure()
    plt.bar(df['disease'], df['size'])
    plt.axhline(y=10, linestyle="--", color="grey", linewidth=2)
    plt.ylabel('LCC Size')
    plt.xticks(rotation=45, ha='right')
    plt.yscale("log")
    savefig('01_lcc_size')
    plt.show()

    # Plot z-score
    plt.figure()
    plt.bar(df['disease'], df['zscore'])
    plt.axhline(y=2, linestyle="--", color="grey", linewidth=2)
    plt.ylabel('LCC Z-Score')
    plt.xticks(rotation=45, ha='right')
    savefig('02_lcc_zscore')
    plt.show()

    # Plot p-value
    plt.figure()
    plt.bar(df['disease'], df['pval'])
    plt.axhline(y=0.05, linestyle="--", color="grey", linewidth=2)
    plt.ylabel('LCC p-value')
    plt.xticks(rotation=45, ha='right')
    plt.yscale("log")
    savefig('03_lcc_pval')
    plt.show()


lcc_size = pd.DataFrame(columns=["disease", "size", "zscore", "pval"])

for d, genes in disease_genes.items():
    data = netmedpy.lcc_significance(ppi, 
                                     genes, 
                                     null_model="log_binning", 
                                     n_iter=5000)
    new_line = [d, data["lcc_size"], data["z_score"], data["p_val"]]
    lcc_size.loc[len(lcc_size.index)] = new_line

save(lcc_size, "output/lcc_size.pkl")
# lcc_size = load("lcc_size.pkl")

plot_lcc_significance(lcc_size)


### Calculate shortest path distance matrix between node pairs
sp_distance = netmedpy.all_pair_distances(ppi, distance="shortest_path", 
                                          n_processors=n_procs, n_tasks=1000)

### CALCULATE PROXIMITY FROM VITAMIN D TO ALL DISEASES
vit_d = {"Vitamin D": targets}
screen_data = netmedpy.screening(vit_d, 
                                 disease_genes, 
                                 ppi, 
                                 sp_distance, 
                                 score="proximity", 
                                 properties=["z_score", "raw_amspl"],
                                 null_model="log_binning", 
                                 n_iter=10000, 
                                 n_procs=n_procs)

save(screen_data, "output/screen.pkl")
# screen_data = load("screen.pkl")


# Results of the screening
def plot_screening(sdata):
    vdata = sdata.T.reset_index()
    vdata.columns = ["disease", "zscore"]
    vdata = vdata.sort_values(by="zscore", ascending=True)

    # Proximity
    plt.figure()
    plt.axhline(y=0, linestyle="--", linewidth=2, color="grey")
    plt.plot(vdata['disease'], vdata['zscore'], marker='o', 
             linestyle='-', markersize=12, linewidth=2)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Proximity Z-Score')
    savefig('04_proximity')
    plt.show()


plot_screening(screen_data['z_score'])


### Closer look at proximity between Vitamin D, Inflammation, and Factor IX deficiency
vdata = screen_data['z_score'].T.reset_index()
vdata.columns = ["disease", "zscore"]
vdata = vdata.sort_values(by="zscore", ascending=True).reset_index(drop=True)

d_name1 = vdata.loc[0, 'disease']
d1 = netmedpy.proximity(ppi, targets, disease_genes[d_name1], sp_distance, 
                        null_model="log_binning", n_iter=10000, symmetric=False)

d_name2 = vdata.loc[len(vdata.index)-1, 'disease']
d2 = netmedpy.proximity(ppi, targets, disease_genes[d_name2], sp_distance, 
                        null_model="log_binning", n_iter=10000, symmetric=False)

save((d_name1, d_name2, d1, d2), "output/d1_d2.pkl")


# Plot results
def plot_histograms(d_name1, d_name2, d1, d2, yl1, yl2, yl3):
    plt.figure(figsize=(10, 7))
    sns.kdeplot(d1['dist'], color='blue', fill=True, 
                alpha=0.5, label=d_name1, bw_adjust=1.5)
    sns.kdeplot(d2['dist'], color='red', fill=True, 
                alpha=0.5, label=d_name2, bw_adjust=1.5)

    plt.axvline(x=d1['raw_amspl'], linewidth=10, color='blue', alpha=0.5)
    plt.text(x=d1['raw_amspl'], y=yl2, s=f"  Z = {d1['z_score']:.2f}", 
             color='blue', verticalalignment='top', horizontalalignment='left', 
             fontsize=20)

    plt.axvline(x=d2['raw_amspl'], linewidth=10, color='red', alpha=0.5)
    plt.text(x=d2['raw_amspl'], y=yl3, s=f"  Z = {d2['z_score']:.2f}", 
             color='red', verticalalignment='top', horizontalalignment='left', 
             fontsize=20)

    plt.xlabel('AMSPL from Vitamin D')
    plt.ylabel('Density')
    plt.ylim(0, yl1)
    plt.legend(loc='upper center', frameon=False)
    savefig('05_amspl_density')
    plt.show()


plot_histograms(d_name1, d_name2, d1, d2, 7, 6.5, 6.5)

### AMSPL equivalent under different distance metrics
amspl = {"Shortest Path": screen_data["raw_amspl"]}

# Random Walks
sp_distance = netmedpy.all_pair_distances(ppi, distance="random_walk")
sd = netmedpy.screening(vit_d, 
                        disease_genes, 
                        ppi, 
                        sp_distance, 
                        score="proximity", 
                        properties=["raw_amspl"], 
                        null_model="log_binning", 
                        n_iter=10, 
                        n_procs=n_procs)

amspl["Random Walks"] = sd["raw_amspl"]

# Biased Random Walks
sp_distance = netmedpy.all_pair_distances(ppi, distance="biased_random_walk")
sd = netmedpy.screening(vit_d, 
                        disease_genes, 
                        ppi, 
                        sp_distance, 
                        score="proximity", 
                        properties=["raw_amspl"], 
                        null_model="log_binning", 
                        n_iter=10, 
                        n_procs=n_procs)

amspl["Biased Random Walks"] = sd["raw_amspl"]

# Communicability
sp_distance = netmedpy.all_pair_distances(ppi, distance="communicability")
sd = netmedpy.screening(vit_d, 
                        disease_genes, 
                        ppi, 
                        sp_distance, 
                        score="proximity", 
                        properties=["raw_amspl"], 
                        null_model="log_binning", 
                        n_iter=10, 
                        n_procs=n_procs)

amspl["Communicability"] = sd["raw_amspl"]

save(amspl, "output/amspl.pkl")
# amspl = load("output/amspl.pkl")

# Plot results
def plot_amspl(amspl, z_score):
    metrics = ["Shortest Path", "Random Walks", "Biased Random Walks", "Communicability"]
    df = amspl["Shortest Path"].T.reset_index()
    df.columns = ["disease", "Shortest Path"]

    for metric in metrics[1:]:
        dnew = amspl[metric].T.reset_index()
        dnew.columns = ["disease", metric]
        df = pd.merge(df, dnew, on="disease")

    df = df.sort_values(by="Shortest Path")

    z_score = screen_data["z_score"]
    zdf = z_score.T.reset_index()
    zdf.columns = ["disease", "z_score"]
    df = pd.merge(df, zdf, on="disease")
    df = df.sort_values(by="z_score")

    plt.figure(figsize=(10, 8))
    for m in metrics:
        plt.plot(df['disease'], df[m] / max(df[m]), marker='o', linestyle='-', 
                 markersize=12, linewidth=2, label=m)

    plt.ylabel("Normalized AMSPL")
    plt.xticks(rotation=90, ha='center')
    plt.legend(loc="best", frameon=False)
    savefig('06_amspl_screening')
    plt.show()

    corr = df.drop('disease', axis=1).astype(float).corr(method='spearman')
    plt.figure(figsize=(10, 8))

    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    sns.heatmap(corr, cmap=cmap, center=0, vmax=1, vmin=-1, annot=True, fmt=".2f", 
                annot_kws={'size': 20, 'color': 'black'}, square=True, 
                linewidths=.5, cbar_kws={"shrink": 1})

    plt.xticks(np.arange(len(corr.columns)) + 0.5, corr.columns, rotation=45, ha='right')
    plt.yticks(np.arange(len(corr.columns)) + 0.5, corr.columns, rotation=0)
    plt.tight_layout()
    savefig('07_metrics_correlation')
    plt.show()


plot_amspl(amspl, screen_data["z_score"])


### NetMedPy Performance evaluation
@ray.remote
def calc_proximity(size, r, mat, ppi):
    c = cronometer()
    nodes = list(ppi.nodes)
    ag_net = p_network.Network(ppi)

    print(f"Size {size} rep {r}")
    a = random.sample(nodes, size)
    b = random.sample(nodes, size)

    c.tick()
    netmedpy.proximity(ppi, a, b, mat, null_model="log_binning", 
                       n_iter=100, bin_size=100)
    c.tock()
    time_new = c.elapsed_seconds

    c.tick()
    ag_net.get_proximity(a, b, bin_size=100, n_iter=100)
    c.tock()
    time_old = c.elapsed_seconds

    print(f"Done Size {size} rep {r}")
    return size, r, time_old, time_new


def size_increasing_time(mat, ppi, ini_size, end_size, increment, num_cpus, reps):
    v = np.arange(ini_size, end_size, increment)
    df = pd.DataFrame({"size": v, "old": 0, "new": 0})
    df.set_index("size", inplace=True)

    ray.shutdown()
    ray.init(num_cpus=num_cpus, log_to_driver=True)

    mat_ref = ray.put(mat)
    ppi_ref = ray.put(ppi)
    futures = []

    for s in v:
        for r in range(reps):
            f = calc_proximity.remote(s, r, mat_ref, ppi_ref)
            futures.append(f)

    res = list(ray.get(futures))
    ray.shutdown()

    for s, r, o, n in res:
        df.loc[s, "old"] += o
        df.loc[s, "new"] += n

    df["old"] = (1.0 * df["old"]) / reps
    df["new"] = (1.0 * df["new"]) / reps
    df.reset_index(inplace=True)
    return df


ppi_file = "data/input/ppi/ppi_network.pkl"
ppi = load(ppi_file)

print("Calculating proximity size dependent")

ini = 50
end = 250
step = 10
reps = 10
num_cpus = n_procs

ptime = size_increasing_time(sp_distance, ppi, ini, end, step, num_cpus, reps)
ptime.to_csv("output/performance_size.csv")

# Plot results
def quadratic(x, a, b, c):
    return a * x ** 2 + b * x + c


params_old, _ = opt.curve_fit(quadratic, ptime["size"], ptime["old"])
ao, bo, co = params_old
x_o = np.linspace(min(ptime["size"]), max(ptime["size"]), 100)
y_o = quadratic(x_o, ao, bo, co)

params_new, _ = opt.curve_fit(quadratic, ptime["size"], ptime["new"])
an, bn, cn = params_new
x_n = np.linspace(min(ptime["size"]), max(ptime["size"]), 100)
y_n = quadratic(x_n, an, bn, cn)

plt.figure(figsize=(10, 7))
plt.scatter(ptime["size"], ptime["old"], color="black", s=60, label="Not optimized")
plt.plot(x_o, y_o, color="black")
plt.scatter(ptime["size"], ptime["new"], color="red", s=60, label="NetMedPy")
plt.plot(x_n, y_n, color="red")

plt.xlabel("Size")
plt.ylabel("log(t) (s)")
plt.yscale("log")
plt.legend(loc="best", frameon=False)
savefig('08_time_execution')
plt.show()
