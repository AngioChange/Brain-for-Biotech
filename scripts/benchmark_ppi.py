import glob
import networkx as nx
import os
import pandas as pd
import re
import shutil
import sys
import useful_functions 

# Reading the training genes
tip_cell_genes = pd.read_csv("G:/Leo/BfBio_reloaded/data/training_genes/Tip_cell_genes_all_tissues_SC_Atlas.csv")
tip_cell_genes = tip_cell_genes["gene"].to_list()

angiogenic_genes = pd.read_csv("G:/Leo/BfBio_reloaded/data/training_genes/AngioGenes_FromSymbols2Entrez.csv")
angiogenic_genes = angiogenic_genes["name"].to_list()

inflammatory_genes = pd.read_csv("G:/Leo/BfBio_reloaded/data/training_genes/training_genes_pro-inflammation.csv")
inflammatory_genes = inflammatory_genes["Gene"].to_list()

brain_genes = pd.read_csv("G:/Leo/BfBio_reloaded/data/training_genes/Markers_brain_vein_2.csv")
brain_genes = brain_genes["gene"].to_list()

lung_genes = pd.read_csv("G:/Leo/BfBio_reloaded/data/training_genes/Markers_lung_vein_1.csv")
lung_genes = lung_genes["gene"].to_list()

# Setting the dataset to process
dataset = sys.argv[1]
if len(sys.argv) != 3:
    sys.exit("Two arguments are required: PPI database and process")

# Setting the process of interest
process = sys.argv[2]
if len(sys.argv) != 3:
    sys.exit("Two arguments are required: PPI database and process")

if process == "inflammation":
    genes = inflammatory_genes
elif process == "angiogenesis":
    genes = angiogenic_genes
elif process == "tip-cell":
    genes = tip_cell_genes
elif process == "lung":
    genes = lung_genes
elif process == "brain":
    genes = brain_genes

# Setting the path to the omics folder 
path_data = f"G:/Leo/BfBio_reloaded/data/annotation-based/{dataset}"

# Create a folder for the upcoming threshold benchmark
path_benchmark = f"G:/Leo/BfBio_reloaded/results/{process}"
path_benchmark_dataset = f"{path_benchmark}/{dataset}"
path_results = f"{path_benchmark_dataset}/Benchmark"

if not os.path.exists(f"{path_benchmark}"):
    os.mkdir(f"{path_benchmark}")

if not os.path.exists(f"{path_benchmark_dataset}"):
    os.mkdir(f"{path_benchmark_dataset}")

if not os.path.exists(f"{path_results}"):
    os.mkdir(f"{path_results}")
else:
    shutil.rmtree(f"{path_results}")
    os.mkdir(f"{path_results}")
   
# Reading the graph file 
G = nx.read_graphml(f"{path_data}/graph_{dataset}.graphml")

# Running a benchmark to determine the best damping factor
benchmark_results = useful_functions.analyse(G, genes)
benchmark_results.to_csv(f"{path_results}/PageRank_evaluation_{dataset}_{process}.csv",
        sep = ",", index = False)

# Retrieving the best DF and saving the file as hyperparameters
hyperparameters_df = benchmark_results[benchmark_results["auc_roc"] == benchmark_results["auc_roc"].max()]
hyperparameters_df.to_csv(f"{path_data}/Hyperparameters_{dataset}_{process}.csv", sep = ",", index = False)
