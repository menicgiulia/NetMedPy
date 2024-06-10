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

#### I. Without installing the package

1. Ensure you have Python installed.
  
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
      
#### II. With installing the package

1. Installing the necessary dependencies:
   
   
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

## Examples

To test the pipeline, we refer to the Vitamin D example, which can be found in the `examples/VitaminD` directory. There are two files that you can use for testing:

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

    
## Package Structure
Root folder organization (__init__.py files removed for simplicity):
```plaintext
│   .gitignore
│   environment.yml                                 // yml file to create conda enviorement
│   README.md
│   setup.py               
│
├───images                                          // directory with figures from paper
│   └───OverviewPipeline.png                        // pipeline flowchart figure from paper
│
└───examples                                        // directory with working examples using the NetMedPy pipeline
│   │   
│   ├───VitamindD                                   // directory with Vitamin D example using the NetMedPy pipeline
│   │    ├───Figures_v2.py                          // python script to recreate the figures from the paper  
│   │    ├───VitD_pipeline.py                       // python script with Vitamin D example using the NetMedPy pipeline  
│   │    ├───VitD_pipeline.ipynb                    // Jupyter notebook with Vitamin D example using the NetMedPy pipeline  
│   │    ├───data                                   // directory with pickle and csv files necessary to get the Vitamin D example working             
│   │    │    ├───Alias.csv                          
│   │    │    ├───disease_genes.pkl                 
│   │    │    ├───ppi_network.pkl                    
│   │    │    └───vitd_targets.pkl
│   │    └───output                                 // directory where the output files from the Vitamin D example are saved
│   │                          
│   └───Basic_example.py                            // python script with dummy data to test the pipeline
│
└───netmedpy                                        // directory containing the python scripts that contain the functions of the NetMedPy pipeline
      ├───DistanceMatrix.py                       
      └───NetMedPy.py
```

## Further information

- Details about each function (what is it used for, what are the input parameters, the possible values of the input parameters, what is the output) from the pipeline are available in the `netmedpy/NetMedPy.py` script in the comments before each function. 
- An example on the use of the implemented functions is available in the file `examples/Basic_example.py', which can be executed fairly quickly in order to test the proper installation of the package and it's functionalities.
- A more elaborate example is available in the files `examples/VitaminD/VitD_pipeline.py` and `examples/VitaminD/VitD_pipeline.ipynb`, testing the functions with different parameters for evaluating the role of Vitamin D in the modulation of
different diseases from a network medicine perspective. The data files (the protein-protein interation network, the disease genes, and the Vitamin D targets) needed for executing this example are available in `examples/VitaminD/data`.

## License

This project is licensed under the terms of the MIT license.


## References

<b id="f1">1</b> Barabási, A. L., Gulbahce, N., & Loscalzo, J. (2011). Network medicine: a network-based approach to human disease. Nature reviews genetics, 12(1), 56-68.[DOI 10.1038/nrg2918](https://doi.org/10.1038/nrg2918) [↩](#a1)

<b id="f2">2</b> Menche, Jörg, et al. "Uncovering disease-disease relationships through the incomplete interactome." Science 347.6224 (2015). [DOI 10.1126/science.1257601](https://doi.org/10.1126/science.1257601) [↩](#a2)

<b id="f3">3</b> Guney, Emre, et al. "Network-based in silico drug efficacy screening." Nature Communications 7,1 (2015). [DOI 10.1038/ncomms10331](https://doi.org/10.1038/ncomms10331) [↩](#a3)
