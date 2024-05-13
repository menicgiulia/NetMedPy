# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 22:07:41 2023

@author: aalda
"""
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import pandas as pd
import random
import time

# Dummy function for demonstration purposes
def evaluate_disease_pair(context):
    # Simulating some processing for this demonstration
    sleep_time = random.uniform(0.1, 1.0)
    time.sleep(sleep_time)
    return ["DiseaseA", "DiseaseB", 0.5, 0.6, 0.1]

def init_pool(CONTEXT):
    global global_context
    global_context = CONTEXT

def worker(k):
    res = evaluate_disease_pair(global_context)
    return res

def eval_k_pairs(k, CONTEXT, outFile):
    global global_context
    global_context = CONTEXT

    with tqdm(total=k, desc="Processing", dynamic_ncols=True, position=0, leave=True) as pbar:
        with Pool(processes=cpu_count(), initializer=init_pool, initargs=(CONTEXT,)) as pool:
            results = []
            for result in pool.imap_unordered(worker, list(range(k))):
                results.append(result)
                pbar.update()

    df = pd.DataFrame(data=results, columns=['D1', 'D2', 'z1', 'z2', 'z1-z2'])
    return df

if __name__ == "__main__":
    CONTEXT = "SomeContext"  # Just a placeholder, you should replace this with your actual context data
    outFile = "results.csv"
    k = 100  # For example, if you want to evaluate 100 pairs

    df = eval_k_pairs(k, CONTEXT, outFile)
