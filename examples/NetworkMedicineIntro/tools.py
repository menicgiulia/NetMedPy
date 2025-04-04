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
