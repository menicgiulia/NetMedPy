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
                        

### Setting up a work environment

#### 1. Without installing the package

1.1. Ensure you have Python installed.
  
1.2. Copy the project at your local or remote machine:
   ```bash
      git clone https://github.com/menicgiulia/NetMedPy.git
   ```
1.3. Navigate to the project directory:
   ```bash
      cd NetMedPy-main
   ```

1.4. Installing the necessary dependencies
   
   
##### Option A: working with Conda

It is recommended to work with Conda, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Create a new conda environment with the environment.yml:

      conda env create -f environment.yml

- Activate your conda :

      conda activate netmedpy_test
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

      pip install networkx seaborn matplotlib numpy pandas ray 
            
1.5. Set up your PYTHONPATH (Replace "/user_path_to/NetMedPy-main/netmedpy" with the appropriate path of the package in your local/remote machine.):

    * On Linux/Mac:
      
      export PYTHONPATH="/user_path_to/NetMedPy-main/netmedpy":$PYTHONPATH
  
    * On Windows shell:
  
      set PYTHONPATH="C:\\user_path_to\\NetMedPy-main\\netmedpy";%PYTHONPATH%

    * On Powershell:

      $env:PYTHONPATH = "C:\\user_path_to\\NetMedPy-main\\netmedpy;" + $env:PYTHONPATH

1.6. Navigate to the directory "examples":
   
      ```bash
      cd examples
      ```
1.7. Run the Basic_example.py script using Python 3 or higher:
   
      ```bash
      python Basic_example.py
      ```
      
#### 2. With installing the package

2.1. Installing the necessary dependencies
   
   
##### Option A: working with Conda

It is recommended to work with Conda, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Download the environment.yml and navigate to the directory of your local/remote machine where the file is located.

- Create a new conda environment with the environment.yml:

      conda env create -f environment.yml

- Activate your conda :

      conda activate netmedpy_test
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

      pip install networkx seaborn matplotlib numpy pandas ray

2.2. Install the package:

      ```bash
      pip install netmedpy
      ```
      
2.3. Download the directory examples.
   
2.4. Navigate to the directory "examples":

      ```bash
      cd /user_path_to/examples
      ```
      
2.5. Run the Basic_example.py script using Python 3 or higher:

      ```bash
      python Basic_example.py
      ```

## Examples

To test the pipeline, we refer to the Vitamin D example, which can be found in the `examples/VitaminD` directory. There are two files that you can use for testing:

- A Python script: `VitD_pipeline.py`
- A Jupyter notebook: `VitD_pipeline.ipynb`

### Instructions

1. Download the `examples` directory:
If you haven't already, download the examples directory from the repository to your local or remote machine. This directory contains all the necessary files to run the example.

2. Prepare the Data:
In the subdirectory VitaminD/data there are the files that contain the necessary data to execute the example, ensure the data files ther. The output files will be stored in the VitaminD/output subdirectory.

3. Navigate to the VitaminD directory:
     ```bash
      cd /user_path_to/examples/VitaminD
      ```
4. Run the Example:

##### Option A: using the Python Script

      python VitD_pipeline.py
      
##### Option B: using the Jupyter Notebook

- Start the Jupyter Kernel

    a) If you are working on a local machine:
    ```bash
      jupyter notebook --browser="browser_of_choise"
    ```
    
  Replace browser_of_choice with your preferred browser (e.g., chrome, firefox). The browser window should pop up automatically. If it doesn't, copy and paste the link provided in the terminal into your browser. The link should look something like this:

    
   * http://localhost:8889/tree?token=5d4ebdddaf6cb1be76fd95c4dde891f24fd941da909129e6
    
       
    b) If you are working on a remote machine:
  
    ```bash
      jupyter notebook --no-browser
    ```
    
  Then copy and paste the link provided in the terminal in your local browser of choise, it should look something like this:

    
   * http://localhost:8888/?token=9feac8ff1d5ba3a86cf8c4309f4988e7db95f42d28fd7772
    
    
- Navigate to the VitD_pipeline.ipynb in the Jupyter Notebook interface and start executing the cells.






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
│   │    ├───Figures_v2.py                          // python script to recreate the figures from the paper  
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



