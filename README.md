# NetMedPy: A Python package for Network Medicine 
#### Authors: Andrés Aldana, Michael Sebek, Gordana Ispirova, Rodrigo Dorantes-Gilardi, Giulia Menichetti (giulia.menichetti@channing.harvard.edu)

## Introduction

Network medicine is a post-genomic discipline that harnesses network science principles to analyze the intricate interactions within biological systems, viewing diseases as localized disruptions in networks of genes, proteins, and other molecular entities <sup id="a1">[1](#f1)</sup>.

The structure of the biological network plays an essential role in the system’s ability to efficiently propagate signals and withstand random failures. Consequently, most analyses in Network Medicine focus on quantifying the efficiency of the communication between different regions of the interactome or protein-protein interaction network.

NetMedPy evaluates network localization (statistical analysis of the largest connected component/subgraph or LCC) <sup id="a2">[2](#f2)</sup>, calculates proximity <sup id="a3">[3](#f3)</sup> and separation <sup id="a2">[2](#f2)</sup> between biological entities, and conducts screenings involving a large number of diseases and drug targets. NetMedPy extends the traditional Network Medicine analyses by providing four default network metrics (shortest paths, random walk, biased random walk, communicability) and four null models (perfect degree match, degree logarithmic binning, strength logarithmic binning, uniform). The user is allowed to introduce custom metrics and null models.


The pipeline workflow is depicted in the figure below.
![Overview Pipeline](https://raw.githubusercontent.com/menicgiulia/NetMedPy/main/images/OverviewPipeline.png)

This Python implementation uses precomputed distance matrices to optimize calculations. With precalculated distances between every node pair, the code can rapidly compute proximity and separation.

## Getting Started      

### Setting up a work environment

NetMedPy has specific requirements for compatibility and ease of use.

### Python Version

NetMedPy requires Python 3.8 or newer, but it is not compatible with Python 3.12 due to incompatibility with Ray. Ensure your Python version is between 3.8 and 3.11.9 inclusive.
    
### Recommended Environment
While not essential, we recommend creating a dedicated conda environment for NetMedPy to ensure all dependencies are properly isolated.

### Required Packages
The following Python packages are required to run NetMedPy:

- Python (>= 3.8, <= 3.11.9)
- numpy
- pandas
- ray
- networkx
- scipy
- matplotlib
- seaborn

### Installation steps

Users can install NetMedPy and its dependencies using PIP (recommended). Alternatively, the source code can be downloaded, allowing for manual installation of the required dependencies if more customization is needed.


#### I. Installing package and dependencies with PIP

##### Option A: working with Conda

Working with Conda is recommended, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Download the environment.yml and navigate to the directory of your local/remote machine where the file is located.

- Create a new conda environment with the `environment.yml` file:

  ```bash
  conda env create -f environment.yml
  ```

- Activate your new conda environment:

  ```bash
  conda activate netmedpy_environment
  ```
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

  ```bash
  pip install networkx seaborn matplotlib numpy pandas ray scipy
  ```

2. Install the package:

  ```bash
  pip install netmedpy
  ```
      
3. Download the directory `examples`.
   
4. Navigate to the directory `examples`:

  ```bash
  cd /user_path_to/examples
  ```
      
5. Run the `Basic_example.py` script using Python 3 or higher:

  ```bash
  python Basic_example.py
  ```
#### II. Manual installation

1. Ensure you have Python >= 3.8, <= 3.11.9 installed.
  
2. Copy the project to your local or remote machine:

  ```bash
  git clone https://github.com/menicgiulia/NetMedPy.git
  ```
3. Navigate to the project directory:

  ```bash
  cd NetMedPy-main
  ```

4. Installing the necessary dependencies:
   
   
##### Option A: working with Conda

Working with Conda is recommended, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Create a new conda environment with the `environment.yml` file:
  
  ```bash
  conda env create -f environment.yml
  ```

- Activate your new conda environment:

  ```bash
  conda activate netmedpy_environment
  ```
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

  ```bash
  pip install networkx seaborn matplotlib numpy pandas ray scipy
  ```
            
5. Set up your PYTHONPATH (Replace `/user_path_to/NetMedPy-main/netmedpy` with the appropriate path of the package in your local/remote machine.):

    _On Linux/Mac_:

    ```bash
    export PYTHONPATH="/user_path_to/NetMedPy-main/netmedpy":$PYTHONPATH
    ```
      
    _On Windows shell_:

    ```bash
    set PYTHONPATH="C:\\user_path_to\\NetMedPy-main\\netmedpy";%PYTHONPATH%
    ```
      
    _On Powershell_:

    ```bash
    $env:PYTHONPATH = "C:\\user_path_to\\NetMedPy-main\\netmedpy;" + $env:PYTHONPATH
    ```
    
6. Navigate to the directory `examples`:
   
  ```bash
  cd examples
  ```
      
7. Run the `Basic_example.py` script using Python 3 or higher (up to 3.11.9, due to conflicts with `Ray`):
   
  ```bash
  python Basic_example.py
  ```

## Example for understanding the pipeline - Vitamin D 

After you have successfully run the Basic_example.py script to further test the pipeline, we refer to the Vitamin D example, which can be found in the `examples/VitaminD` directory. There are two files that you can use for testing:

- A Python script: `VitD_pipeline.py`
- A Jupyter notebook: `VitD_pipeline.ipynb`

### Instructions on testing the Vitamin D example

1. Download the `examples` directory:
If you haven't already done so, download the `examples` directory from the repository to your local or remote machine. This directory contains all the necessary files to run the examples.

2. Prepare the Data:
In the subdirectory `VitaminD/data` there are the files that contain the necessary data to execute the example, ensure the data files there. The output files will be stored in the `VitaminD/output` subdirectory.

3. Navigate to the `VitaminD` directory:
 
  ```bash
  cd /user_path_to/examples/VitaminD
  ```
     
4. Run the Example:

##### Option A: using the Python Script

  ```bash
  python VitD_pipeline.py
  ```
      
##### Option B: using the Jupyter Notebook

- Make sure you have the `jupyter` package installed.

  ```bash
  pip install jupyter
  ```
  
- Start the Jupyter Kernel

    a) If you are working on a local machine:
  
    ```bash
    jupyter notebook --browser="browser_of_choice"
    ```
    
  Replace `browser_of_choice` with your preferred browser (e.g., chrome, firefox). The browser window should pop up automatically. If it doesn't, copy and paste the link provided in the terminal into your browser. The link should look something like this:

    
   * http://localhost:8889/tree?token=5d4ebdddaf6cb1be76fd95c4dde891f24fd941da909129e6
    
       
    b) If you are working on a remote machine:
  
    ```bash
    jupyter notebook --no-browser
    ```
    
  Then copy and paste the link provided in the terminal in your local browser of choice. It should look something like this:

    
   * http://localhost:8888/?token=9feac8ff1d5ba3a86cf8c4309f4988e7db95f42d28fd7772
    
    
- Navigate to the `VitD_pipeline.ipynb` in the Jupyter Notebook interface and start executing the cells.


### Extract and evaluate disease modules

  - From a dictionary of diseases `disease_genes` the function lcc_significance will calculate the statistical significance of the size of the Largest Connected Component (LCC) of a subgraph induced by the node set `genes` in the network `ppi`. This function generates a null model distribution for the LCC size by resampling nodes from the network while preserving their degrees (`null_model="log_binning"`). The statistical significance of the observed LCC size is then determined by comparing it against this null model distribution.
  
  - The parameter `null_model` can be `degree_match`, `log_binning`, `uniform`, or `custom` (defined by the user).
  
```python
#Load disease genes dictonary from the pickle file in `examples/VitaminD/data/disease_genes.pkl`
with open("examples/VitaminD/data/disease_genes.pkl","rb") as file:
  disease_genes = pickle.load(file)

lcc_size = pd.DataFrame(columns = ["disease","size","zscore","pval"])

for d,genes in disease_genes.items():
    data = netmedpy.lcc_significance(ppi, genes,
                                     null_model="log_binning",n_iter=10000)
    new_line = [d,data["lcc_size"],data["z_score"],data["p_val"]]
    lcc_size.loc[len(lcc_size.index)] = new_line

#Keep only diseases with an LCC larger than 10 and statistically significant
#Filtering the disease sets to the LCC is optional and not mandatory for the subsequent analyses
significant = lcc_size.query("size > 10 and zscore > 2 and pval<0.05")
disease_names = significant.disease
```
### Evaluate Average Minimum Shortest Path Length (AMSPL) between Vitamin D and Inflammation and between Vitamin D and Factor IX Deficiency disease

  - The function proximity calculates the proximity between two sets of nodes in a given graph based on the approach described by Guney et al., 2016. The method computes either the average minimum shortest path length (AMSPL) or its symmetrical version (SASPL) between two sets of nodes.

   In this example, the function calculates the proximity between the Vitamin D targets stored in `examples/VitaminD/data/vitd_targets.pkl` and the disease genes from the `examples/VitaminD/data/disease_genes.pkl` file for the two diseases: `Inflammation` and `Factor IX Deficiency`. The null model of choice, in this case, is `log_binning`.

   - The function returns a dictionary containing various statistics related to proximity, including:
       - 'd_mu': The average distance in the randomized samples.
       - 'd_sigma': The standard deviation of distances in the randomized samples.
       - 'z_score': The z-score of the actual distance in relation to the randomized samples.
       - 'p_value_single_tail': One-tail P-value associated with the proximity z-score
       - 'p_value_double_tail': Two-tail P-value associated with the proximity z-score
       - 'p_val': P-value associated with the z-score.
       - 'raw_amspl': The raw average minimum shortest path length between the two sets of interest.
       - 'dist': A list containing distances from each randomization iteration.

       
```python
#Load PPI network
with open("examples/VitaminD/data/ppi_network.pkl","rb") as file:
  ppi = pickle.load(file)

#Load drug targets
with open("examples/VitaminD/data/vitd_targets.pkl","rb") as file:
  targets = pickle.load(file)

#Load disease genes
with open("examples/VitaminD/data/disease_genes.pkl","rb") as file:
  disease_genes = pickle.load(file)

inflammation = netmedpy.proximity(ppi, targets,
                                  dgenes["Inflammation"], sp_distance,
                                  null_model="log_binning",n_iter=10000,
                                  symmetric=False)

factorix = netmedpy.proximity(ppi, targets,
                                  dgenes["Factor IX Deficiency"], sp_distance,
                                  null_model="log_binning",n_iter=10000,
                                  symmetric=False)

plot_histograms(inflammation, factorix)
```

### Evaluate Average Minimum Shortest Path Length (AMSPL) under different distances 

- The function `all_pair_distances` calculates distances between every pair of nodes in a graph according to the specified method and returns a DistanceMatrix object. This function supports multiple distance calculation methods, including shortest path, various types of random walks, and user-defined methods.
  
- The function `screening` screens for relationships between sets of source and target nodes within a given network, evaluating proximity or separation. This function facilitates drug repurposing and other network medicine applications by allowing the assessment of network-based relationships.

- In this example using the `all_pair_distances` function the distance between every pair of nodes in the protein-protein interaction network stored in the file `examples/VitaminD/data/ppi_network.pkl` are calculated, using different parameters for the method of calculation: `random_walk`, `biased_random_walk`, and `communicability`.

- For each calculation of the distance matrix the AMSPL is calculated using the `screening` function evaluating `proximity`.

```python
#Load PPI network
with open("examples/VitaminD/data/ppi_network.pkl","rb") as file:
  ppi = pickle.load(file)

#Load drug targets
with open("examples/VitaminD/data/vitd_targets.pkl","rb") as file:
  targets = pickle.load(file)

#Load disease genes
with open("examples/VitaminD/data/disease_genes.pkl","rb") as file:
  disease_genes = pickle.load(file)

#Shortest Paths
amspl = {"Shortest Path":screen_data["raw_amspl"]}

#Random Walks
sp_distance = netmedpy.all_pair_distances(ppi,distance="random_walk")
screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                 sp_distance,score="proximity",
                                 properties=["raw_amspl"],
                                 null_model="log_binning",
                                 n_iter=10,n_procs=20)

amspl["Random Walks"] = screen_data["raw_amspl"]

#Biased Random Walks
sp_distance = netmedpy.all_pair_distances(ppi,distance="biased_random_walk")
screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                 sp_distance,score="proximity",
                                 properties=["raw_amspl"],
                                 null_model="log_binning",
                                 n_iter=10,n_procs=20)

amspl["Biased Random Walks"] = screen_data["raw_amspl"]


#Communicability
sp_distance = netmedpy.all_pair_distances(ppi,distance="communicability")
screen_data = netmedpy.screening(vit_d, dgenes, ppi,
                                 sp_distance,score="proximity",
                                 properties=["raw_amspl"],
                                 null_model="log_binning",
                                 n_iter=10,n_procs=20)

amspl["Communicability"] = screen_data["raw_amspl"]
```
## Example for entry-level users - Introduction to Network Medicine

This example introduces the core concepts of network medicine through a guided analysis of Vitamin D's relationship to several diseases using protein-protein interaction networks. The Jupyter notebook (`Intro_Network_Medicine.ipynb`) provides a step-by-step workflow demonstrating how to build and analyze biological networks to uncover drug-disease relationships.

1. Clone or update the repository

```bash
git clone https://github.com/menicgiulia/NetMedPy.git    # if you haven't already
cd NetMedPy
git pull                                                # if you already cloned
```

2. Inspect the input data
   
Pre-downloaded disease gene lists are located at:

 ```plaintext
examples/NetworkMedicineIntro/input_data/disease_genes/
├── DGN_Huntington.csv      # Huntington's disease gene list
├── DGN_Rickets.csv         # Rickets gene list
├── DGN_VDdeff.csv          # Vitamin D deficiency gene list
└── DGN_inflammation.csv    # Inflammation gene list
```

(You do not need to manually download STRING data—this is handled by the notebook.)

3. Install prerequisites
   
Activate your NetMedPy conda environment or install via pip:

```bash
# If using conda:
conda activate netmedpy_environment

# Ensure NetMedPy and Jupyter are installed:
pip install netmedpy jupyter
```

4. Navigate to the example folder
   
```bash
cd examples/NetworkMedicineIntro
```

5. Launch and run the notebook
   
```bash
jupyter notebook Intro_Network_Medicine.ipynb
```

Execute cells in order—each cell saves outputs under examples/NetworkMedicineIntro/output/.

### Notebook workflow - Steps in `Intro_Network_Medicine.ipynb`:

**1. Download and filter STRING PPI data** 

The notebook first defines the URL for STRING v12 and downloads the protein-protein interaction data:

```python
# Define the URL for the STRING PPI dataset
string_url = "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz"

# Define paths for temporary files
string_gz_path = './tmp_string/string.gz'

# Download and extract STRING data
print("Downloading STRING dataset...")
tools.download_file(string_url, string_gz_path)
tools.ungz_file(string_gz_path, "./tmp_string/string_data")
```

It then processes the data by removing prefixes and converting Ensembl IDs to HGNC symbols:

```python
print("Processing protein names...")
string_df["protein1"] = string_df["protein1"].str.replace("9606.", "", regex=False)
string_df["protein2"] = string_df["protein2"].str.replace("9606.", "", regex=False)

# Convert Ensembl IDs to HGNC symbols
ens_to_hgnc = tools.ensembl_to_hgnc(string_df)
string_df["HGNC1"] = string_df["protein1"].map(ens_to_hgnc)
string_df["HGNC2"] = string_df["protein2"].map(ens_to_hgnc)
```

Finally, it filters the network, extracts the largest connected component, and saves it:

```python
filtered_df = string_df.query("weight > 300")
G_string = nx.from_pandas_edgelist(filtered_df, 'HGNC1', 'HGNC2', create_using=nx.Graph())

G_string = netmedpy.extract_lcc(G_string.nodes, G_string)

# Save to CSV
df_edges = nx.to_pandas_edgelist(G_string)
df_edges.to_csv("output/string_ppi_filtered.csv", index=False)
```



**2. Extract Vitamin D targets**
   
The notebook extracts compound-protein databases from a pre-packaged zip file:

```python
# Define database directory path
data_path = "./output/cpie_Databases"

if os.path.exists(data_path):
    shutil.rmtree(data_path)

tools.unzip_file("../VitaminD/supplementary/sup_data/cpie_databases/Databases.zip", data_path)
```

It then loads multiple databases into memory and searches for Vitamin D (Cholecalciferol) targets:

```python
# Store all databases in a dictionary
dbs = {
    'chembl': chembl_data,
    'bdb': BDB_data,
    'stitch': sttch_data,
    'ctd': CTD_data,
    'dtc': DTC_data,
    'db': DB_data,
    'dc': DC_data
}

# Cholecalciferol (PubChem CID: 5280795)
comp_id = 5280795

# Initialize Comp2Prot
C2P = Comp2Prot('local', dbs=dbs)

# Search for interactions
comp_dat, status = C2P.comp_interactions(input_id=comp_id)

# Extract HGNC symbols
vd_targets = {"Vitamin D": list(comp_dat.hgnc_symbol)} 

# Save extracted targets
with open('./output/vd_targets.json', 'w') as f:
    json.dump(vd_targets, f)
```

**3. Load disease gene sets**

Disease-gene associations are loaded from DisGeNet files and filtered by confidence score:

```python
# Directory containing the disease genes
dis_gene_path = "input_data/disease_genes"

disease_file_names = {
    "Huntington":"DGN_Huntington.csv",
    "Inflammation": "DGN_inflammation.csv",
    "Rickets": "DGN_Rickets.csv",
    "Vit. D defficiency": "DGN_VDdeff.csv"
}

disease_genes = {}

# Load files and filter for strong associations
for name,file_name in disease_file_names.items():
    path = dis_gene_path + "/" + file_name

    df = pd.read_csv(path)
    df = df.query("Score_gda > 0.1")

    disease_genes[name] =  list(df.Gene)

# Save file
with open('./output/disease_genes.json', 'w') as f:
    json.dump(disease_genes, f)
```

**4.  Verify network coverage**
   
The notebook checks which disease genes and drug targets are found in the PPI network:

```python
# Load PPI network
ppi = pd.read_csv("output/string_ppi_filtered.csv")
ppi = nx.from_pandas_edgelist(ppi, 'source', 'target', create_using=nx.Graph())

# Keep only associations existing in the PPI
nodes = set(ppi.nodes)
for name, genes in disease_genes.items():
    disease_genes[name] = set(genes) & nodes
    print(f"{name}: {len(disease_genes[name])} associations in PPI")

for name, targets in dtargets.items():
    dtargets[name] = set(targets) & nodes
    print(f"{name}: {len(dtargets[name])} targets in PPI")
```

**5. Compute random walk distances**
    
The notebook calculates biased random walk distances between all nodes:


```python
# Calculate Random Walk based distance between all pair of genes
dmat = netmedpy.all_pair_distances(
    ppi,
    distance='biased_random_walk',
    reset = 0.3
)

# Save distances for further use
netmedpy.save_distances(dmat,"output/ppi_distances_BRW.pkl")
```

**6. Calculate proximity with log-binning null model**

The notebook computes proximity z-scores using the log-binning null model:

```python
# Calculate proximity between Vitamin D targets and Diseases
proximity_lb = netmedpy.screening(
    dtargets, 
    disease_genes, 
    ppi,
    dmat,
    score="proximity",
    properties=["z_score"],
    null_model="log_binning",
    n_iter=10000,n_procs=10
)

zscore_lb = proximity_lb['z_score'].T
zscore_lb = zscore_lb.sort_values(by='Vitamin D')
zscore_lb
```

**7. Repeat analysis with degree-matched null model**

The same analysis is performed using the degree-match null model for comparison:

```python
pythonproximity_dm = netmedpy.screening(
    dtargets, 
    disease_genes, 
    ppi,
    dmat,
    score="proximity",
    properties=["z_score"],
    null_model="degree_match",
    n_iter=10000,n_procs=10
)

zscore_dm = proximity_dm['z_score'].T
zscore_dm = zscore_dm.sort_values(by='Vitamin D')
zscore_dm
```


**8. Compare results from both null models**

Finally, the notebook combines results from both methods for comparison:
pythonzscore_lb.columns = ["Log Binning"]
zscore_dm.columns = ["Degree Match"]

zscore = pd.merge(zscore_lb,zscore_dm, left_index=True, right_index=True)

zscore
This produces a table showing z-scores from both null models, with Vitamin D deficiency having the strongest connection to Vitamin D targets.



### Data sources

- STRING v12: Human protein-protein interactions downloaded directly from stringdb-downloads.org
- Compound-target databases: Collection of databases accessed from the VitaminD supplementary data folder
- DisGeNet: Disease-gene associations provided as CSV files in the input_data folder

### Expected outputs

After running the notebook, the following files will be created:


```bash
output/
├── string_ppi_filtered.csv    # Filtered STRING PPI network
├── vd_targets.json           # Vitamin D protein targets
├── disease_genes.json        # Disease gene sets
├── ppi_distances_BRW.pkl     # Biased random walk distance matrix
└── cpie_Databases/           # Extracted compound-protein interaction databases
```

## Package Structure
Root folder organization (__init__.py files removed for simplicity):
```plaintext
│   .gitattributes                                 
│   .gitignore                                    
│   LICENSE.txt                                     // License information for the package
│   README.md                                       // Package documentation
│   environment.yml                                 // yml file to create conda environment
│   setup.py                                        // Package installation script
│
├───doc                                             // Documentation directory
│   └───source                                      // Source files for documentation
│       │   DistanceMatrix.rst                      // Documentation for DistanceMatrix module
│       │   NetMedPy.rst                            // Documentation for NetMedPy module
│       │   conf.py                                 // Sphinx configuration file for documentation
│       │   index.rst                               // Main index file for documentation
│       │   Makefile                                // Make file for building documentation
│       │   make.bat                                // Batch script for building documentation on Windows
│
├───examples                                        // directory with working examples using the NetMedPy pipeline
│   │   Basic_example.py                            // python script for running a basic example to test the pipeline
│   │   Cronometer.py                               // Performance timing utility
│   │   VitD_pipeline.ipynb                         // Jupyter notebook with Vitamin D example using the NetMedPy pipeline
│   │   VitD_pipeline.py                            // python script with Vitamin D example using the NetMedPy pipeline
│   │   1_4_netsize_edges.png                       // Figure showing network size and edges relationships
│   │   1_7_prox_vd.png                             // Figure related to proximity and Vitamin D
│   │   1_8_correlation.png                         // Correlation analysis figure
│   │   2_2_deviation.png                           // Deviation analysis figure
│   │   2_3_rank_correlation_distr...               // Rank correlation distribution figure
│   │
│   ├───NetworkMedicineIntro                        // Introduction to Network Medicine examples
│   │   │   Intro_Network_Medicine.ipynb            // Jupyter notebook with intro to network medicine
│   │   │   tools.py                                // Helper tools for the analysis
│   │   │
│   │   └───input_data/disease_genes                // Disease gene data for examples
│   │           DGN_Huntington.csv                  // Huntington disease gene data
│   │           DGN_Rickets.csv                     // Rickets disease gene data
│   │           DGN_VDdeff.csv                      // Vitamin D deficiency gene data
│   │           DGN_inflammation.csv                // Inflammation gene data
│   │
│   └───VitaminD                                    // directory with Vitamin D example using the NetMedPy pipeline
│       ├───data                                    // directory with data files necessary for the Vitamin D example
│       │   └───input                               // Input data directory
│       │       ├───disease_genes                   // Disease gene data directory
│       │       │       disease_genes_merge.pkl     // Merged disease genes data
│       │       │
│       │       ├───drug_targets                    // Drug target data directory
│       │       │       vitd_targets_cpie.pkl       // Vitamin D targets data
│       │       │
│       │       └───ppi                             // Protein-protein interaction data
│       │               ppi_network.pkl             // PPI network data
│       │               Alias.csv                   // Alias mapping file
│       │
│       ├───guney                                   // Implementation of Guney's network algorithms
│       │       distances.py                        // Distance calculation functions
│       │       network.py                          // Network manipulation functions
│       │
│       ├───output                                  // directory where the output files from the Vitamin D example are saved
│       │       amspl.pkl                           // Analysis output file
│       │       d1_d2.pkl                           // Disease pairs data
│       │       inf_fix.pkl                         // Inflammation-related output
│       │       inf_hun.pkl                         // Huntington-related output
│       │       lcc_size.pkl                        // Largest connected component size data
│       │       performance_size.csv                // Performance metrics
│       │       screen.pkl                          // Screening results
│       │
│       └───supplementary                           // Supplementary materials
│           └───sup_code                            // Supplementary code
│               └───data_integration                // Data integration scripts
│
├───images                                          // directory with figures from paper
│       OverviewPipeline.png                        // pipeline flowchart figure from paper
│
└───netmedpy                                        // directory containing the python scripts that contain the functions of the NetMedPy pipeline
        DistanceMatrix.py                           // Module for distance matrix calculations
        NetMedPy.py                                 // Core NetMedPy functionality
```

## Further information

- Details about each function (what is it used for, what are the input parameters, the possible values of the input parameters, what is the output) from the pipeline are available in the `netmedpy/NetMedPy.py` script in the comments before each function. 
- An example on the use of the implemented functions is available in the file `examples/Basic_example.py', which can be executed fairly quickly in order to test the proper installation of the package and its functionalities.
- A more elaborate example is available in the files `examples/VitaminD/VitD_pipeline.py` and `examples/VitaminD/VitD_pipeline.ipynb`, testing the functions with different parameters for evaluating the role of Vitamin D in the modulation of
different diseases from a network medicine perspective. The data files (the protein-protein interaction network, the disease genes, and the Vitamin D targets) needed for executing this example are available in `examples/VitaminD/data`.

## License

This project is licensed under the terms of the MIT license.


## References

<b id="f1">1</b> Barabási, A. L., Gulbahce, N., & Loscalzo, J. (2011). Network medicine: a network-based approach to human disease. Nature Reviews Genetics, 12(1), 56-68.[DOI 10.1038/nrg2918](https://doi.org/10.1038/nrg2918) [↩](#a1)

<b id="f2">2</b> Menche, Jörg, et al. "Uncovering disease-disease relationships through the incomplete interactome." Science 347.6224 (2015). [DOI 10.1126/science.1257601](https://doi.org/10.1126/science.1257601) [↩](#a2)

<b id="f3">3</b> Guney, Emre, et al. "Network-based in silico drug efficacy screening." Nature Communications 7,1 (2015). [DOI 10.1038/ncomms10331](https://doi.org/10.1038/ncomms10331) [↩](#a3)
