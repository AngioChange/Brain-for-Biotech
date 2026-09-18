import glob
import networkx as nx
import os
import pandas as pd
import re
import shutil
import sys
import omics_analysis
from tqdm import tqdm

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
    sys.exit("Two arguments are required: dataset and process")

# Setting the process of interest
process = sys.argv[2]
if len(sys.argv) != 3:
    sys.exit("Two arguments are required: dataset and process")

if process != "inflammation":
    process_dir = "angio_tip-cell"
else:
    process_dir = "inflammation" 

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
path_data = f"G:/Leo/BfBio_reloaded/data/graphs/omics/{process_dir}"

# Setting the path to the hyperparameters folder
path_hyperparameters = f"G:/Leo/BfBio_reloaded/data/graphs/omics/hyperparameters/{process}/{dataset}"
if not os.path.exists(f"{path_hyperparameters}"):
    os.mkdir(f"{path_hyperparameters}")

# Listing all the graph files to consider for the input dataset 
files = glob.glob(f"{path_data}/{dataset}*.graphml")

# Create a folder for the upcoming threshold benchmark
path_benchmark = f"G:/Leo/BfBio_reloaded/results/{process}"
path_benchmark_dataset = f"{path_benchmark}/{dataset}"
path_results = f"{path_benchmark_dataset}/Benchmark_correlation_and_DF"

if not os.path.exists(f"{path_benchmark}"):
    os.mkdir(f"{path_benchmark}")

if not os.path.exists(f"{path_benchmark_dataset}"):
    os.mkdir(f"{path_benchmark_dataset}")

if not os.path.exists(f"{path_results}"):
    os.mkdir(f"{path_results}")
else:
    shutil.rmtree(f"{path_results}")
    os.mkdir(f"{path_results}")
   
# Initializing the workflow    
for file in tqdm(files, desc = "Processing files ..."):

# Reading the network file and retrieving the correlation threshold that has been applied to generate the graph file
    G = nx.read_graphml(file)
    filename = os.path.basename(file)
    match = re.search(f"{dataset}_(.*)_graph.graphml", filename)

    if match:
        correlation_threshold = match.group(1)
        correlation_threshold = correlation_threshold.replace(".", "")

# Determining the output file name
    output_file = f"{path_results}/benchmark_{dataset}_corr_{correlation_threshold}.csv"

# Running a benchmark on the damping factor value to find the best
    benchmark_results = omics_analysis.analyse(G, genes, correlation_threshold, output_file)    

# Analyzing the benchmark results : we need to determine the best threshold for the scores and the best damping factor
correlation_thresholds = []
AUCs = []
Damping_factors = []

files = glob.glob(f"{path_results}/*.csv")
for file in files:
    results = pd.read_csv(file)
    correlation_thresholds.append(int(results.loc[:, "correlation_threshold"].mean()))

    best_df = results[results["auc_roc"] == results["auc_roc"].max()]
    for i in best_df.itertuples():
        Damping_factors.append(i[1])
        AUCs.append(i[7])

benchmark_process = pd.DataFrame({"correlation_threshold": correlation_thresholds,
    "Best_AUC": AUCs,
    "Best_DF": Damping_factors})

benchmark_process.to_csv(f"{path_results}/Benchmark_damping_factor_{dataset}_{process}.csv",
    sep = ",", index = False)

# Using the benchmark to find the hyperparameters (correlation threshold and damping factor)
hyperparameters_df = benchmark_process[benchmark_process["Best_AUC"] == benchmark_process["Best_AUC"].max()]
hyperparameters_df.to_csv(f"{path_hyperparameters}/Hyperparameters_{dataset}_{process}.csv", sep = ",", index = False)
