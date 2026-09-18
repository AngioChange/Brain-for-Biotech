import os
import useful_functions
from datetime import datetime
import networkx as nx
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

# String dataset
info_file = "9606.protein.info.v12.0.txt"
link_file = "9606.protein.links.v12.0.txt"
string_dir = "G:/Leo/BfBio_reloaded/data/annotation-based/STRING"
training_genes_dir = "G:/Leo/BfBio_reloaded/data/training_genes"

def create_out_dir():
    today_date = datetime.today().strftime('%Y_%m_%d')  # Format: YYYY-MM-DD
    out_dir = f"/data/Julia/String_output_{today_date}/"
    # Create the folder if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Folder '{out_dir}' created successfully!")
    else:
        print(f"Folder '{out_dir}' already exists.")
    return out_dir
    
def make_network():
    string_network = pd.read_csv(f"{string_dir}/{link_file}" , sep = " ")
    protein_info = pd.read_csv(f"{string_dir}/{info_file}" , sep = "\t")
    id_to_name = dict(zip(protein_info["#string_protein_id"], protein_info["preferred_name"]))
    # Add gene names to the network data
    string_network["gene1"] = string_network["protein1"].map(id_to_name)
    string_network["gene2"] = string_network["protein2"].map(id_to_name)
    string_network = string_network[string_network['gene1'] != string_network['gene2']]
    # Remove rows where mapping is not possible (e.g., missing gene names)
    string_network = string_network.dropna(subset=["gene1", "gene2"])
    return string_network

def make_string_graph(string_network , threshold):
    filtered_data = string_network[string_network["combined_score"] > threshold]
    max_score = filtered_data["combined_score"].max()
    filtered_data["normalized_score"] = filtered_data["combined_score"] / max_score
    print(f"shape of filtered data: {filtered_data.shape}")

    string_graph = nx.Graph()  # Directed graph for PageRank
    edges = filtered_data[["gene1", "gene2", "normalized_score"]].values
    string_graph.add_weighted_edges_from(edges)
    return string_graph

def analyse(genes, out_dir, process):
    string_network = make_network()
    thresholds = [0, 100, 200, 300, 400, 500, 600, 700]
    for threshold in thresholds:
        print(f"Processing graph for threshold{threshold}")
        string_graph = make_string_graph(string_network , threshold)
        n_nodes = len(string_graph.nodes())
        n_edges = len(string_graph.edges())
        valid_genes = [node for node in string_graph.nodes(data=False) if node in genes]  # Extract valid genes

        kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Ensure reproducibility
        folds = list(kf.split(valid_genes))  # Get indices for splits

        all_results = []
        k = 100
        for df in np.round(np.arange(0.1,1, 0.1), 1):
            fold_evaluations = useful_functions.run_page_rank_databases(df, string_graph, folds, valid_genes , k)  # Function 
    
            # Collect results for each fold
            for eval_result in fold_evaluations:
                # Add damping factor and fold to the result
                eval_result['damping_factor'] = df
                # Append to the results list
                all_results.append(eval_result)


 
        results_df = pd.DataFrame(all_results)
        average_results_string_per_damping_factor = results_df.groupby("damping_factor").mean()
        average_results_string_per_damping_factor['n_nodes'] = n_nodes
        average_results_string_per_damping_factor['n_edges'] = n_edges
        average_results_string_per_damping_factor['n_tr_genes'] = len(valid_genes)
        average_results_string_per_damping_factor['threshold'] = threshold
        # Reset index for cleaner display
        average_results_string_per_damping_factor = average_results_string_per_damping_factor.reset_index()
    
        
        output_file = f"{out_dir}/page_rank_evaluations_string_graph_{threshold}_{process}.csv"
        average_results_string_per_damping_factor.to_csv(output_file, index=False)
    
    
def main():
    """
    Main function of the script.
    """
 
    tip_cell_genes=pd.read_csv(f"{training_genes_dir}Tip_cell_genes_all_tissues_SC_Atlas.csv")
    tip_cell_genes = tip_cell_genes['gene'].to_list()
    print(f"number of tip cell genes {len(tip_cell_genes)}" )
    angiogenic_genes = pd.read_csv(f"{training_genes_dir}AngioGenes_FromSymbols2Entrez.csv")
    angiogenic_genes = angiogenic_genes['name'].to_list()
    print(f"number of angiogenic genes {len(angiogenic_genes)}" )
    
    print(f"running analysis for tip cell genes")
    analyse(tip_cell_genes)
    
    print(f"running analysis for angiogenic genes")
    analyse(angiogenic_genes)

if __name__ == "__main__":
    main()