# NetMedPy: A Python package for Network Medicine 
#### Authors: Andrés Aldana, Michael Sebek, Gordana Ispirova, Rodrigo Dorantes-Gilardi, Giulia Menichetti (giulia.menichetti@channing.harvard.edu)

## Introduction

Network medicine is a post-genomic discipline that harnesses network science principles to analyze the intricate interactions within biological systems, viewing diseases as localized disruptions in networks of genes, proteins, and other molecular entities [1].

The structure of the biological network plays an essential role in the system’s ability to efficiently propagate signals and withstand random failures. Consequently, most analyses in Network Medicine focus on quantifying the efficiency of the communication between different regions of the interactome or protein-protein interaction network.

NetMedPy evaluates network localization (statistical analysis of the largest connected component/subgraph or LCC [2]), calculates proximity [3] and separation [2] between
biological entities, and conducts screenings involving a large number of diseases and drug targets. NetMedPy extends the traditional Network Medicine analyses by providing four default network metrics (shortest paths, random walk, biased random walk, communicability) and four null models (perfect degree match, degree logarithmic binning, strength logarithmic binning, uniform). The user is allowed to introduce custom metrics and null models.

The pipeline workflow is depicted in the figure below.
![Pipeline](/images/OverviewPipeline.png)

This Python implementation uses precomputed distance matrices to optimize calculations. With precalculated distances between every node pair, the code can rapidly compute proximity and separation.

## Getting Started
      
                            
```
### Setting up a work environment

## Without installing the package

- Ensure you have Python installed.
- Copy the project at your local or remote machine:

git clone https://github.com/menicgiulia/NetMedPy.git

- Navigate to the project directory:

cd NetMedPy-main

It is recommended to work with Conda, but it is not essential. If you chose to work with Conda, these are the steps you need to take:
- Ensure you have Conda installed.
- Create a new conda environment with the environment.yml:

conda env create -f environment.yml

- Activate your conda :

conda activate netmedpy_test

- Set up your PYTHONPATH:

  On Linux/Mac:
      
      export PYTHONPATH="/user_path_to/NetMedPy/netmedpy":$PYTHONPATH
  
  On Windows shell:
  
      set PYTHONPATH="C:\\user_path_to\\NetMedPy\\netmedpy";%PYTHONPATH%

  On Powershell:

      $env:PYTHONPATH = "C:\\user_path_to\\NetMedPy\\netmedpy;" + $env:PYTHONPATH


- Install necessary libraries:

  ```bash
  pip install networkx seaborn matplotlib scipy
  ```

Import NetworkMetrics module and networkx:

```python
import random
import networkx
import network_distances.NetworkMetrics as distances
```

Replace with the appropriate path to the `NetworkMetrics` module on your system.

### Create the network and sets

As an example, consider a basic network G, and two sets of randomly selected node sets:

```python
G = nx.barabasi_albert_graph(1000, m)
nodes = list(G.nodes())
S = random.sample(nodes, 20)
T = random.sample(nodes, 20)
```

### Precomputing distances

```python
 D = distances.all_pair_distances(G)
```

### Save distance matrix as a pickle file

```python
 distances.save_distances(D, "distances.pkl")
```

### Load distance matrix from a pickle file

```python
D = distances.load_distances("distances.pkl")
```

### Calculate proximity between S and T, using exact degree preserving randomization and 1000 iterations

```python
p = distances.proximity(G, S, T, D, degree_preserving='exact', n_iter=1000)
```

The result is a dictionary with the next attributes:

- 'd_mu': The average distance in the randomized samples.
- 'd_sigma': The standard deviation of distances in the randomized samples.
- 'z_score': The z-score of the actual distance in relation to the randomized samples.
- 'raw_amspl': The raw average minimum shortest path length between sets T and S.
- 'dist': A list containing distances from each randomization iteration.

To retrieve the z-score value of proximity:

```python
print(p['z_score'])
```

### Calculate separation between S and T

Calculate the separation between S and T:

```python
sep = distances.separation(G, S, T, D)
print(sep)
```

### Calculate the z-score significance of the separation value:

```python
s = distances.separation_z_score(G, S, T, D, degree_preserving='exact', n_iter=1000)
print(s['z_score'])
```
### Calculate LCC significance

```python
#Let's select a random set of nodes to evaluate:
G1 = G.subgraph(random.sample(nodes, 50))
L = set(G1.nodes())

lcc_data = metrics.lcc_significance(G,L, degree_preserving='exact',n_iter=10000)

# Z-score and p-value
size = lcc_data['lcc_size']
z = lcc_data['z_score']
p = lcc_data['p_val']

print(f"LCC-size={size} z-score={z:0.2f} p-value={p:0.2f}")
```

## Package Structure
Root folder organization (__init__.py files removed for simplicity):
```plaintext
│   .gitignore
│   environment.yml                                 // yml file to create conda enviorement
│   README.md
│
├───images                                          // directory with figures from paper
│   └───OverviewPipeline.png                        // pipeline flowchart figure from paper
│
└───examples                                        // directory with working examples using the netmedpy pipeline
│   │   
│   ├───VitamindD                                   // directory with Vitamin D example using the netmedpy pipeline      
│   │    ├───VitD_pipeline.py                       // python script with Vitamin D example using the netmedpy pipeline  
│   │    ├───VitD_pipeline.ipynb                    // Jupyter notebook with Vitamin D example using the netmedpy pipeline  
│   │    ├───data                                   // directory with pickle and csv files necessary to get the Vitamin D example working             
│   │    │    ├───Alias.csv                          
│   │    │    ├───disease_genes.pkl                 
│   │    │    ├───ppi_network.pkl                    
│   │    │    └───vitd_targets.pkl
│   │    └───output                                 // directory where the output files from the Vitamin D example are saved
│   │                          
│   └───Basic_example.py                            // python script with dummy data to test the pipeline
│
└───netmedpy                                        // directory containing the python scripts that contain the functions of the netmedpy pipeline
      ├───DistanceMatrix.py                       
      └───NetMedPy.py
```

## Further reading

An example on the use of the implemented functions is available in the file 'Example.py'. Consult the full documentation of the appropiate functions in the file 'NetworkMetrics.pdf' or 'NetworkMetrics.md'

## References

[1] Barabási, A. L., Gulbahce, N., & Loscalzo, J. (2011). Network medicine: a network-based approach to human disease. Nature reviews genetics, 12(1), 56-68.[DOI 10.1038/nrg2918](https://doi.org/10.1038/nrg2918)
[2] Menche, Jörg, et al. "Uncovering disease-disease relationships through the incomplete interactome." Science 347.6224 (2015). [DOI 10.1126/science.1257601](https://doi.org/10.1126/science.1257601)
[3] Guney, Emre, et al. "Network-based in silico drug efficacy screening." Nature Communications 7,1 (2015). [DOI 10.1038/ncomms10331](https://doi.org/10.1038/ncomms10331)



