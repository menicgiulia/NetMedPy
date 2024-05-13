# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 19:54:26 2024

@author: aalda
"""

import multiprocessing
import ray
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define the computation task
def computation_task(x):
    # Example task: compute the square and simulate a workload
    return np.random.normal(loc=0, scale=5, size=1000000)

# Function to measure multiprocessing performance
def measure_multiprocessing_performance(tasks):
    start_time = time.time()
    with multiprocessing.Pool() as pool:
        pool.map(computation_task, range(tasks))
    return time.time() - start_time

# Function to measure Ray performance
@ray.remote
def ray_computation_task(x):
    # This is essentially the same as computation_task but adapted for Ray
    return np.random.normal(loc=0, scale=5, size=10000)

def measure_ray_performance(tasks):
    ray.init(ignore_reinit_error=True)
    start_time = time.time()
    ray.get([ray_computation_task.remote(i) for i in range(tasks)])
    ray.shutdown()
    return time.time() - start_time

# Main comparison function
def compare_frameworks(task_counts):
    mp = []
    ra = []

    for tasks in task_counts:
        mp_time = measure_multiprocessing_performance(tasks)
        ray_time = measure_ray_performance(tasks)

        mp.append(mp_time)
        ra.append(ray_time)

    return pd.DataFrame({"Tasks":task_counts,"mp_time":mp,"ray_time":ra})

# Example usage
if __name__ == "__main__":
    task_counts = [100,200,300,400,500,600,700,800,900,1000]  # Adjust based on your testing needs
    df = compare_frameworks(task_counts)


    plt.figure(figsize=(9,6))
    plt.plot(df.Tasks,df.mp_time,label="MP",lw=2,color="black")
    plt.plot(df.Tasks,df.ray_time,label="Ray",lw=2,color="red")
    plt.xlabel("Tasks",fontsize=14)
    plt.ylabel("Time (s)",fontsize=14)
    plt.legend(loc="best")
    plt.show()
