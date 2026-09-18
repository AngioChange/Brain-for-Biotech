import os
import networkx as nx
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
import useful_functions
from datetime import datetime

annot_graph_dir = "G:/Leo/BfBio_reloaded/data/graphs/annotation-based/BioGrid/"
training_genes_dir = "G:/Leo/BfBio_reloaded/data/training_genes/"
annot_network_dir = "G:/Leo/Thèse/BfBio_reloaded/data/annotation-based/BioGrid/"

def create_out_dir():
    today_date = datetime.today().strftime('%Y_%m_%d')  # Format: YYYY-MM-DD
    out_dir = f"C:/Users/betto/Desktop/Léo/Thèse/BfBio_reloaded/data/Annotation-based/BioGrid/BioGrid_output_{today_date}/"
    # Create the folder if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Folder '{out_dir}' created successfully!")
    else:
        print(f"Folder '{out_dir}' already exists.")
    return out_dir

def get_largest_cc(genes):
    biogrid_filtered = pd.read_csv(f"{annot_network_dir}biogrid_network.csv")

    # Step 2: Create an undirected graph
    BioGridGraph = nx.Graph()  # Use Graph() for an undirected graph

    # Step 3: Add edges to the graph
    # Iterate through the rows of the DataFrame and add edges with weights
    for index, row in biogrid_filtered.iterrows():
        BioGridGraph.add_edge(row['source'], row['target'], weight=row['confidence'])
    
    print(f"The graph is directed? {BioGridGraph.is_directed()}")
    if nx.is_connected(BioGridGraph):
        print("The graph is fully connected!")
    else:
        print("The graph is not fully connected!")
    
    # Find the connected components and their sizes
    components = list(nx.connected_components(BioGridGraph))
    component_sizes = [len(component) for component in components]
    max_component_size = max(component_sizes)
    # Print the number of connected components and their sizes
    print(f"The graph has {len(components)} connected components.")
    print(f"max component size is {max_component_size}")
    
    non_connected_components = [component for component in components if len(component) < max_component_size]

    # Flatten the list of nodes
    nodes_in_non_connected_components = set().union(*non_connected_components)
    print(f"Nodes in non-connected components: {len(nodes_in_non_connected_components)}")
    
    genes = set(genes)
    overlap_genes = nodes_in_non_connected_components.intersection(genes)
    print(f"Genes in non-connected components that are also in the training set: {len(overlap_genes)}")
    print("These genes are:", overlap_genes)
    
    largest_cc = max(components, key=len)
    largest_cc_subgraph = BioGridGraph.subgraph(largest_cc)

    # Verify the size of the largest connected component
    print(f"Largest connected component has {len(largest_cc_subgraph.nodes)} nodes and {len(largest_cc_subgraph.edges)} edges.")
    return largest_cc_subgraph

def analyse(genes, out_dir, process):
    biogrid_graph = get_largest_cc(genes)
    n_nodes = len(biogrid_graph.nodes())
    n_edges = len(biogrid_graph.edges())
    valid_genes = [node for node in biogrid_graph.nodes(data=False) if node in genes]  # Extract valid genes
  
   
    # Step 2: Split genes into 5 folds
    kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Ensure reproducibility
    folds = list(kf.split(valid_genes))  # Get indices for splits
    
    all_results = []
    k = 100
    for df in np.round(np.arange(0.1,1, 0.1), 1):
        fold_evaluations = useful_functions.run_page_rank_databases(df, biogrid_graph, folds , valid_genes , k ) 
    
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

    # Display the averaged results
    if process == "angiogenesis":
        file_name = f"{out_dir}/page_rank_evaluations_biogrid_graph_angio.csv"
    elif process == "tip-cell":
        file_name = f"{out_dir}/page_rank_evaluations_biogrid_graph_tip_cell.csv"
    elif process == "inflammation":
        file_name = f"{out_dir}/page_rank_evaluations_biogrid_graph_inflammation.csv"
    
    average_results_per_damping_factor.to_csv(file_name, index = False)

    
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
    
    out_dir = create_out_dir()
    analyse(tip_cell_genes , out_dir)
    analyse(angiogenic_genes , out_dir)
    
if __name__ == "__main__":
    main()


def run_page_rank_biogrid(df, graph , k, folds, valid_genes):
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
