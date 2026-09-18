import string_analysis

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

# Setting the process of interest
process = sys.argv[1]
if len(sys.argv) != 2:
    sys.exit("The process was not specified !")

# Setting the path to the STRING data
path_data = "G:/Leo/BfBio_reloaded/data/annotation-based/STRING"

# Setting the path to the ConsensusPath graph
path_graph = "G:/Leo/BfBio_reloaded/data/graphs/annotation-based/STRING"


# Create a folder for the upcoming threshold benchmark
path_benchmark = f"G:/Leo/BfBio_reloaded/results/{process}"
path_benchmark_string = f"{path_benchmark}/STRING"
path_results = f"{path_benchmark_string}/Benchmark_thresholds"

if not os.path.exists(f"{path_benchmark}"):
    os.mkdir(f"{path_benchmark}")

if not os.path.exists(f"{path_benchmark_string}"):
    os.mkdir(f"{path_benchmark_string}")

if not os.path.exists(f"{path_results}"):
    os.mkdir(f"{path_results}")
   
    

# Reading the network file
string_network = pd.read_csv(f"{path_data}/string_network.csv")
       
# Running a benchmark on the combined_score values to filter the data
# and find the best threshold for gene prioritization
thresholds = [100, 200, 300, 400, 500, 600, 700]
for threshold in tqdm(thresholds, position = 0, desc = "Processing thresholds ..."):
    print(f"Processing graph for threshold {threshold} and process {process}")

    # Selecting the training genes that will be used in the PageRank algorithm
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


    # Initiating the workflow
    string_graph = string_analysis.make_string_graph(string_network, threshold)
    n_nodes = len(string_graph.nodes())
    n_edges = len(string_graph.edges())
    valid_genes = [node for node in string_graph.nodes(data = False) if node in genes]

    kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
    folds = list(kf.split(valid_genes))

    all_results = []
    k = 100
    for df in np.arange(0.1, 1, 0.1).round(1):
        print("########## Processing damping factor ", df, "##########")
        fold_evaluations = useful_functions.run_page_rank_databases(df, 
            string_graph, 
            folds, 
            valid_genes, 
            k)

        # Retrieving the results for each fold
        for eval_result in fold_evaluations:
            # Add damping factor and fold to the results
            eval_result["damping_factor"] = df
            # Append this to the results list
            all_results.append(eval_result)

    results_df = pd.DataFrame(all_results)
    average_results_string_per_damping_factor = results_df.groupby("damping_factor").mean()
    average_results_string_per_damping_factor["n_nodes"] = n_nodes
    average_results_string_per_damping_factor["n_edges"] = n_edges
    average_results_string_per_damping_factor["n_tr_genes"] = len(valid_genes)
    average_results_string_per_damping_factor["threshold"] = threshold

    # Reset index for a cleaner display
    average_results_string_per_damping_factor = average_results_string_per_damping_factor.reset_index()

    # Saving the results
    average_results_string_per_damping_factor.to_csv(f"{path_results}/pagerank_evaluations_string_graph_{threshold}_{process}.csv",
                                                            sep = ",", index = False)

# Analyzing the benchmark results : we need to determine the best threshold for the scores and the best damping factor
thresholds = []
AUCs = []
Damping_factors = []

files = glob.glob(f"{path_results}/*.csv")
for file in files:
    results = pd.read_csv(file)
    thresholds.append(int(results.loc[:, "threshold"].mean()))

    best_df = results[results["auc_roc"] == results["auc_roc"].max()]
    for i in best_df.itertuples():
        Damping_factors.append(i[1])
        AUCs.append(i[7])

benchmark_process = pd.DataFrame({"Threshold": thresholds,
    "Best_AUC": AUCs,
    "Best_DF": Damping_factors})

benchmark_process.to_csv(f"{path_results}/Benchmark_thresholds_STRING_{process}.csv",
    sep = ",", index = False)

# Using the benchmark to find the hyperparameters (scores threshold and damping factor)
hyperparameters_df = benchmark_process[benchmark_process["Best_AUC"] == benchmark_process["Best_AUC"].max()]
hyperparameters_df.to_csv(f"{path_data}/Hyperparameters_STRING_{process}.csv", sep = ",", index = False)

# Using the best threshold to create a graph file
for i in hyperparameters_df.itertuples():
    best_threshold = i[1]

string_graph = string_analysis.make_string_graph(string_network, threshold)
nx.write_graphml(string_graph, f"{path_graph}/graph_STRING_{process}.graphml")