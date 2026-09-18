import os
import useful_functions
from datetime import datetime
from itertools import combinations
import networkx as nx
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold


training_genes_dir = "G:/Leo/BfBio_reloaded/data/training_genes"
consensus_network_file = "G:/Leo/BfBio_reloaded/data/annotation-based/ConsensusPath/ConsensusPath_network.csv"

def create_out_dir():
    today_date = datetime.today().strftime('%Y_%m_%d')  # Format: YYYY-MM-DD
    out_dir = f"G:/Leo/BfBio_reloaded/data/annotation-based/ConsensusPath/ConsensusPath_output_{today_date}/"
    # Create the folder if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Folder '{out_dir}' created successfully!")
    else:
        print(f"Folder '{out_dir}' already exists.")
    return out_dir

def create_graph(consensus_network_file):
    network = pd.read_csv(consensus_network_file)
    max_score = network["weight"].max()
    network = network[network['source'] != network['target']]

    graph = nx.Graph()  # undirected graph for PageRank
    edges = network[["source", "target", "weight"]].values
    graph.add_weighted_edges_from(edges)
    return graph

def analyse(genes, graph, out_dir):
    n_nodes = len(graph.nodes())
    n_edges = len(graph.edges())
    valid_genes = [node for node in graph.nodes(data=False) if node in genes]  # Extract valid genes

    # Step 2: Split genes into 5 folds
    kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Ensure reproducibility
    folds = list(kf.split(valid_genes))  # Get indices for splits
    len(valid_genes)
    
    all_results = []
    k = 100
    for df in np.round(np.arange(0.1,1, 0.1), 1):
        fold_evaluations = useful_functions.run_page_rank_databases(df, graph, folds, valid_genes , k )  
    
        # Collect results for each fold
        for eval_result in fold_evaluations:
            # Add damping factor and fold to the result
            eval_result['damping_factor'] = df
            # Append to the results list
            all_results.append(eval_result)


    # Convert the results into a Pandas DataFrame
    results_df = pd.DataFrame(all_results)
    average_results_per_damping_factor = results_df.groupby("damping_factor").mean()
    average_results_per_damping_factor['n_nodes'] = n_nodes
    average_results_per_damping_factor['n_edges'] = n_edges
    average_results_per_damping_factor['n_tr_genes'] = len(valid_genes)
    # Reset index for cleaner display
    average_results_per_damping_factor = average_results_per_damping_factor.reset_index()
    if len(valid_genes) <= 47:
        output_file = f"{out_dir}/page_rank_evaluations_consensus_graph_angio.csv"
    else:
        output_file = f"{out_dir}/page_rank_evaluations_consensus_graph_tip_cell.csv"
        
    
    average_results_per_damping_factor.to_csv(output_file, index=False)
    
def main():
    """
    Main function of the script.
    """
    tip_cell_genes=pd.read_csv(f"{training_genes_dir}/Tip_cell_genes_all_tissues_SC_Atlas.csv")
    tip_cell_genes = tip_cell_genes['gene'].to_list()
    print(f"number of tip cell genes {len(tip_cell_genes)}" )
    
    angiogenic_genes = pd.read_csv(f"{training_genes_dir}/AngioGenes_FromSymbols2Entrez.csv")
    angiogenic_genes = angiogenic_genes['name'].to_list()
    print(f"number of angiogenic genes {len(angiogenic_genes)}" )
    graph = create_graph()
    print(f"The GeneCard graph is directed? {graph.is_directed()}")
    if nx.is_connected(graph):
        print("The GeneCard graph is fully connected!")
    else:
        print("The GeneCard graph is not fully connected!")
        
    out_dir = create_out_dir()
    print(f"running analysis for tip cell genes")
    analyse(tip_cell_genes, graph, out_dir)
    
    print(f"running analysis for angiogenic genes")
    analyse(angiogenic_genes, graph, out_dir)

if __name__ == "__main__":
    main()


def create_binary_interactions(participants):
    """
    Function to create binary interactions
    """
    interactors = participants.split(",")
    return list(combinations(interactors, 2))


def run_page_rank_consensus(df, graph , k, folds, valid_genes):
    fold_evaluations = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        fold = fold_idx + 1
        print(f"Processing Fold {fold}...")
        train_genes = [valid_genes[i] for i in train_idx]  # Training set (seed nodes)
        test_genes = [valid_genes[i] for i in test_idx]    # Test set
    
        # Build the personalization vector (1 for training genes, 0 for others)
        personalization = {node: 0 for node in graph.nodes}
        for node  in graph.nodes(data=False):
            if node in train_genes:
                personalization[node] = 1.0
    
        # Normalize the personalization vector
        total = sum(personalization.values())
        personalization = {k: v / total for k, v in personalization.items()}
    
        # Run Personalized PageRank
        ppr_scores = nx.pagerank(graph, alpha=df, personalization=personalization , weight="weight")
        
        ppr_scores_filtered = {
            node: score
            for node, score in ppr_scores.items()
            if node not in train_genes # Ensure the node not in the training genes
        }
        evaluation = useful_functions.evaluate_fold(test_genes, ppr_scores_filtered, k=k)
        fold_evaluations.append({"fold": fold, **evaluation})
    return fold_evaluations