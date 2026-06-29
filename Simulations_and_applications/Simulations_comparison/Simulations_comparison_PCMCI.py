# Simulations for the state-of-the-art PCMCI methods. This code uses the Tigramite package
# Performance will be compared with our method, which is implemented in R. The code for our method can be found in the corrsponding R file.

# generate_data is the function that generated the data with sample size n, number of variables p (p=m from the manuscript) and structure='VAR' or 'GARCH' and heavy_tailed=TRUE or FALSE. 
# For generating random graph we use the function nx.erdos_renyi_graph(p, 1/p, directed=True) from igraph library
# We generate data with a random graph + estimates the graph using one of the methods + compute the distance between true graph and estimated graph + repeat 100 times and return the mean of the distances
# Results are saved for later processing.

import os, sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import time
import networkx as nx
# from concurrent.futures import ProcessPoolExecutor


import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from scipy.stats import pareto, cauchy, norm


import tigramite as tg
from tigramite.independence_tests.robust_parcorr import RobustParCorr
from tigramite.pcmci import PCMCI

import tigramite
from tigramite import data_processing as pp
from tigramite.toymodels import structural_causal_processes as toys

from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.lpcmci import LPCMCI

from tigramite.independence_tests.parcorr import ParCorr
from tigramite.independence_tests.robust_parcorr import RobustParCorr
from tigramite.independence_tests.parcorr_wls import ParCorrWLS
from tigramite.independence_tests.gpdc import GPDC
from tigramite.independence_tests.cmiknn import CMIknn
from tigramite.independence_tests.cmisymb import CMIsymb
from tigramite.independence_tests.gsquared import Gsquared
from tigramite.independence_tests.regressionCI import RegressionCI


## ==== PARAMETERS ====

# Scenario parameters
structure = 'VAR' # 'VAR' or 'GARCH'
heavy_tailed = True # True or False


# Define the ranges for n and p
n_values_cor = [500, 5000] # [500, 5000]
p_values_cor = [2, 4, 7, 10, 20] # [2, 4, 7, 10, 20]
nb_sim_cor = 100 # 100

n_values_gpdc = [500] # [500]
p_values_gpdc = [2, 4, 7] # [2, 4, 7]
nb_sim_gpdc = 100 # 100

# Folder where results will be stored
results_folder = './Results/Simulations_comparison/tables/' # folder where results will be stored

os.makedirs(results_folder, exist_ok=True)

np.random.seed(0)


def generate_data(n=200, p=4, structure="VAR", heavy_tailed=True, seed=None):

    def generate_data_internal(n, p, structure, heavy_tailed, seed):
        # Directed Erdos-Renyi graph, no self-loops
        true_graph = nx.erdos_renyi_graph(p, 1/p, seed=seed, directed=True)
        true_graph.remove_edges_from(nx.selfloop_edges(true_graph))
        
        X = np.zeros((n, p))

        if structure == "VAR":
            ar_coef = 0.1
            cross_coef = 0.5
            noise_sd = 1
            normalize_by_indegree = False

            if heavy_tailed:
                epsilon = pareto.rvs(1, size=(n, p))
            else:
                epsilon = norm.rvs(size=(n, p))

            # Create B[i, k] = effect of source k on target i
            B = np.zeros((p, p))
            for source, target in true_graph.edges():
                B[target, source] = cross_coef

            if normalize_by_indegree:
                indeg = np.sum(B != 0, axis=1)
                indeg[indeg == 0] = 1
                B = B / indeg[:, None]

            np.fill_diagonal(B, ar_coef)
            
            for ti in range(1, n):
                X[ti, :] = B @ X[ti - 1, :] + epsilon[ti, :]

            X = pd.DataFrame(X)
        
        
        
        elif structure == "GARCH":
            if heavy_tailed:
                epsilon = cauchy.rvs(size=(n, p))
            else:
                epsilon = norm.rvs(size=(n, p))

            effect = 0.5

            # A[i, k] = 1 if edge from k to i exists
            A = np.zeros((p, p))
            for source, target in true_graph.edges():
                A[target, source] = 1
            
            for ti in range(1, n):
                addition = (A @ (X[ti - 1, :] ** 2)) * effect
                X[ti, :] = np.sqrt(0.1 + addition) * epsilon[ti, :]

            X = pd.DataFrame(X)
        
        else:
            raise ValueError("structure must be either 'VAR' or 'GARCH'")

        return {"data": X, "true_graph": true_graph}
    
    tmp_seed = np.random.randint(int(1e6))
    X = generate_data_internal(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed, seed=tmp_seed)
    trynm = 0
    while X["data"].iloc[-1, :].isna().sum() != 0:
        tmp_seed = np.random.randint(int(1e6))
        X = generate_data_internal(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed, seed=tmp_seed)
        # DEBUG CHECK: to be sure it is numerically stationary.
        trynm += 1
        if trynm > 100:
            msg = f"Non-stationary DGP (n = {n}, p = {p}) -- Failed to generate valid data after 100 attempts."
            print(msg)
            # warnings.warn(msg)
            raise RuntimeError(msg)

    return {"data": X["data"], "true_graph": X["true_graph"]}






def get_edges_from_graph(graph):
    max_node = max(graph.keys()) + 1
    adj_matrix = [[0] * max_node for _ in range(max_node)]

    # Fill the adjacency matrix based on the edges in the graph
    for node, edges in graph.items():
        for edge in edges:
            dest_node = edge[0]
            if node != dest_node:  # Exclude loop edges
                adj_matrix[node][dest_node] = 1

    edges = []
    # Iterate over the adjacency matrix to find edges
    for i in range(len(adj_matrix)):
        for j in range(len(adj_matrix[i])):
            if adj_matrix[i][j] == 1:
                edges.append((j, i))

    return edges


def distance_between_two_graphs(graph1, graph2):
    l = len(list(set(graph1.edges).intersection(graph2.edges)))
    return len(graph1.edges) + len(graph2.edges) - 2 * l


# Define the function for plotting two directed graphs side by side
def plot_directed_graphs_side_by_side(graph1, graph2):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    nx.draw(graph1, ax=axes[0], with_labels=True, node_color='lightblue', node_size=500, arrows=True)
    axes[0].set_title('Graph 1')

    nx.draw(graph2, ax=axes[1], with_labels=True, node_color='lightgreen', node_size=500, arrows=True)
    axes[1].set_title('Graph 2')

    plt.tight_layout()
    plt.show()







def PCMCI_RobustParCorr_one_simulation(n=1000, p=4,  structure = 'VAR',  heavy_tailed = False):
    data = generate_data(n, p, structure, heavy_tailed)
   
    true_graph = data['true_graph']
    column_names = [f'X{i}' for i in range(1, data['data'].shape[1] + 1)]
    data = pd.DataFrame(data['data'].values, columns=column_names)
   
    tau_max = 1
    pcmci = PCMCI(dataframe=pp.DataFrame(data.values), cond_ind_test=RobustParCorr())
    results = pcmci.run_pcmci(tau_max=tau_max, tau_min=1, pc_alpha=None)
    graph = pcmci.all_parents
   
    estimated = nx.DiGraph(get_edges_from_graph(graph))
    true_graph = nx.DiGraph(list(true_graph.edges))
   
    return distance_between_two_graphs(estimated, true_graph)


def PCMCI_GPDC_one_simulation(n=1000, p=4,  structure = 'VAR',  heavy_tailed = False):
    data = generate_data(n, p, structure, heavy_tailed)
   
    true_graph = data['true_graph']
    column_names = [f'X{i}' for i in range(1, data['data'].shape[1] + 1)]
    data = pd.DataFrame(data['data'].values, columns=column_names)
   
    tau_max = 1
    pcmci = PCMCI(dataframe=pp.DataFrame(data.values), cond_ind_test=GPDC())
    results = pcmci.run_pcmci(tau_max=tau_max, tau_min=1, pc_alpha=None)
    graph = pcmci.all_parents
   
    estimated = nx.DiGraph(get_edges_from_graph(graph))
    true_graph = nx.DiGraph(list(true_graph.edges))
   
    return distance_between_two_graphs(estimated, true_graph)



script_start_time = time.perf_counter()

if(heavy_tailed):
    ht_str = 'heavy'
else:
    ht_str = 'light'



print("Starting simulations for PCMCI with RobustParCorr...")

# Initialize a list to store results
results_PCMCI_Cor = []

# Loop over all combinations of n and p and store every repetition's output
for n in n_values_cor:
    print(f"-- Starting simulations for n={n} (at {(time.perf_counter() - script_start_time)/60:.2f} minutes).")
    for p in p_values_cor:
        print(f"---- Starting simulations for p={p} (at {(time.perf_counter() - script_start_time)/60:.2f} minutes).")
        for rep in range(1, nb_sim_cor+1):
            out = PCMCI_RobustParCorr_one_simulation(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed)
            results_PCMCI_Cor.append({'n': n, 'p': p, 'replication': rep, 'distance': out})



# Convert results to a DataFrame
results_PCMCI_Cor = pd.DataFrame.from_dict(results_PCMCI_Cor)
results_PCMCI_Cor['method'] = 'PCMCI_Cor'
results_PCMCI_Cor['structure'] = structure
results_PCMCI_Cor['heavy_tailed'] = heavy_tailed
results_PCMCI_Cor = results_PCMCI_Cor[['method', 'structure', 'heavy_tailed', 'n', 'p', 'replication', 'distance']]
results_PCMCI_Cor.to_csv(f'{results_folder}results_PCMCI_Cor_{structure}_{ht_str}.csv', index=False)
results_PCMCI_Cor


del results_PCMCI_Cor




print("Starting simulations for PCMCI with GPDC...")

# Initialize a list to store results
results_PCMCI_GPDC = []

# Loop over all combinations of n and p and store every repetition's output
for n in n_values_gpdc:
    print(f"-- Starting simulations for n={n} (at {(time.perf_counter() - script_start_time)/60:.2f} minutes).")
    for p in p_values_gpdc:
        print(f"---- Starting simulations for p={p} (at {(time.perf_counter() - script_start_time)/60:.2f} minutes).")
        for rep in range(1, nb_sim_gpdc+1):
            out = PCMCI_GPDC_one_simulation(n=n, p=p, structure=structure, heavy_tailed=heavy_tailed)
            results_PCMCI_GPDC.append({'n': n, 'p': p, 'replication': rep, 'distance': out})



# Convert results to a DataFrame
results_PCMCI_GPDC = pd.DataFrame.from_dict(results_PCMCI_GPDC)
results_PCMCI_GPDC['method'] = 'PCMCI_GPDC'
results_PCMCI_GPDC['structure'] = structure
results_PCMCI_GPDC['heavy_tailed'] = heavy_tailed
results_PCMCI_GPDC = results_PCMCI_GPDC[['method', 'structure', 'heavy_tailed', 'n', 'p', 'replication', 'distance']]
results_PCMCI_GPDC.to_csv(f'{results_folder}results_PCMCI_GPDC_{structure}_{ht_str}.csv', index=False)
results_PCMCI_GPDC

del results_PCMCI_GPDC




script_end_time = time.perf_counter()
tot_exec_time = script_end_time - script_start_time
print(f"Script execution time: {tot_exec_time/3600:.2f} hours ({tot_exec_time/60:.2f} minutes).")
# write time to a file
with open(f'{results_folder}execution_time_PCMCI_{structure}_{ht_str}.txt', 'w') as f:
    f.write(f"Script execution time: {tot_exec_time/3600:.2f} hours ({tot_exec_time/60:.2f} minutes). \n")






# ######Time measurement - how long does it take to estimate the graph? 
# data = generate_data(n=500, p=20, structure = 'VAR', heavy_tailed=True)
# column_names = [f'X{i}' for i in range(1, data['data'].shape[1] + 1)]
# data = pd.DataFrame(data['data'].values, columns=column_names)

# start_time = time.time()

# pcmci = PCMCI(dataframe=pp.DataFrame(data.values), cond_ind_test=GPDC()) #Change cond_ind_test=RobustParCorr() for other measurment
# results = pcmci.run_pcmci(tau_max=1, tau_min=1, pc_alpha=None)
# graph = pcmci.all_parents
    
# end_time = time.time()

# end_time - start_time


