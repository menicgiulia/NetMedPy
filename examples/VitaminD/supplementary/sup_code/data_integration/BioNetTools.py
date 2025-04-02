# %%
import pandas as pd
import numpy as np
import requests
import zipfile
import gzip
import networkx as nx
import shutil
import mygene
import os

def ensembl_to_hgnc(df):
    # Initialize Server and data
    mg = mygene.MyGeneInfo()
    unique_ensembl_ids = pd.concat([df["protein1"], df["protein2"]]).unique()

    # Query mygene to map Ensembl IDs to HGNC symbols
    results = mg.querymany(unique_ensembl_ids, scopes="ensembl.protein", fields="symbol", species="human")

    # Create a dictionary mapping Ensembl IDs to HGNC symbols
    ensembl_to_hgnc = {item["query"]: item.get("symbol", "Unknown") for item in results}

    return ensembl_to_hgnc


def download_file(url, save_path):
    try:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        # Send a GET request to the URL
        response = requests.get(url, stream=True)
        # Check if the request was successful
        response.raise_for_status()
        
        # Open the file in binary write mode
        with open(save_path, 'wb') as file:
            # Write the content in chunks to avoid memory issues
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
                
        print(f"File downloaded successfully and saved to {save_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download the file: {e}")
    except OSError as e:
        print(f"Failed to create directory or write file: {e}")

def unzip_file(zip_path, extract_to):
    """
    Unzip a zip file to the specified directory.
    """
    try:
        # Ensure the extraction directory exists
        os.makedirs(extract_to, exist_ok=True)
        
        # Extract the zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        print(f"Files extracted to: {extract_to}")
    except zipfile.BadZipFile as e:
        print(f"Failed to extract the zip file: {e}")


def ungz_file(gz_path, extract_to):
    """
    Extract a gz file to the specified directory.
    """
    try:
        # Ensure the extraction directory exists
        os.makedirs(extract_to, exist_ok=True)
        
        # Determine the output file name by stripping the .gz extension
        base_name = os.path.basename(gz_path)
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]
        output_path = os.path.join(extract_to, base_name)
        
        # Open the gz file and write the decompressed data to the output file
        with gzip.open(gz_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())
                
        print(f"File extracted to: {output_path}")
    except OSError as e:
        print(f"Failed to extract the gz file: {e}")


def load_mitab_to_dataframe(file_path):
    """
    Load a MITAB file (tab-separated) into a pandas DataFrame.
    """
    try:
        # MITAB files are tab-separated; we assume the first row is the header
        df = pd.read_csv(file_path, sep="\t", low_memory=False)
        print(f"MITAB file loaded into DataFrame: {file_path}")
        return df
    except Exception as e:
        print(f"Failed to load MITAB file into DataFrame: {e}")
        return None

def extract_hgnc_biogrid(alt_ids):
    for part in alt_ids.split('|'):
        if 'entrez gene/locuslink:' in part:
            return part.split(':')[-1].strip()  # Extract gene name after the last colon
    return None  # Return None if no gene name is found


def extract_score_biogrid(value):
    if value == '-':
        return 0  
    elif 'score:' in value:
        return float(value.split(':')[1])  # Extract the numeric value and convert to float
    else:
        print(f"Unexpected format {value}")
        return 0  # Default to 0 if unexpected format

def download_string_ppi(url='default', file_path = None):
    url_d = ""
    if url=='default':
        url_d = "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz"
    else:
        url_d = url

    if file_path == None:
        path = './tmp_string/string.gz'
        download_file(url_d,path)
    else:
        path = file_path

    ungz_file(path, "./tmp_string/string_data")

    ppi_df = pd.read_csv("./tmp_string/string_data/string",sep="\s")

    shutil.rmtree("./tmp_string")

    # remove prefix from protein names
    ppi_df["protein1"] = ppi_df["protein1"].str.replace("9606.", "", regex=False)
    ppi_df["protein2"] = ppi_df["protein2"].str.replace("9606.", "", regex=False)

    ens_to_hgnc = ensembl_to_hgnc(ppi_df)

    # Add HGNC1 and HGNC2 columns using the mapping
    ppi_df["HGNC1"] = ppi_df["protein1"].map(ens_to_hgnc)
    ppi_df["HGNC2"] = ppi_df["protein2"].map(ens_to_hgnc)

    # Remove unknowns
    ppi_df = ppi_df.query("HGNC1 != 'Unknown'")
    ppi_df = ppi_df.query("HGNC2 != 'Unknown'")

    return ppi_df

def download_biogrid_ppi(url='default', file_path = None):
    url_d = ""
    if url=='default':
        url_d = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.241/BIOGRID-ALL-4.4.241.mitab.zip"
    else:
        url_d = url

    if file_path == None:
        path = './tmp_biogrid/biogrid.zip'
        download_file(url_d,path)
    else:
        path = file_path

    unzip_file(path, "./tmp_biogrid/biogrid_data")

    df = pd.read_csv("./tmp_biogrid/biogrid_data/BIOGRID-ALL-4.4.241.mitab.txt", sep="\t", low_memory=False)

    shutil.rmtree("./tmp_biogrid")

    df = df.query("`Taxid Interactor A` == 'taxid:9606'")
    df = df.query("`Taxid Interactor B` == 'taxid:9606'")
    df = df[df["Interaction Types"].str.contains("physical association",case=False)]

    df = df[['Alt IDs Interactor A','Alt IDs Interactor B', 
            'Aliases Interactor A', 'Aliases Interactor B','Confidence Values']]

    df['Gene_Name_A'] = df['Alt IDs Interactor A'].apply(extract_hgnc_biogrid)
    df['Gene_Name_B'] = df['Alt IDs Interactor B'].apply(extract_hgnc_biogrid)
    df = df.query("Gene_Name_A != Gene_Name_B")    

    df['Score'] = df['Confidence Values'].copy().apply(extract_score_biogrid)

    df = df[['Gene_Name_A','Gene_Name_B','Score']]

    return df
