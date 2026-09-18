import os
import networkx as nx
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
import useful_functions
from datetime import datetime

omics_graph_dir = "G:/Leo/BfBio_reloaded/data/graphs/omics/"
training_genes_dir = "G:/Leo/BfBio_reloaded/data/training_genes/"

def analyse(graph, genes, correlation_threshold, out_file):
    if nx.is_connected(graph):
        print("The graph is fully connected!")
    else:
        print("The graph is not fully connected.")
        graph = exctract_largest_cc(graph, genes)
    n_nodes = len(graph.nodes())
    n_edges = len(graph.edges())
    valid_genes = [data['name'] for _, data in graph.nodes(data=True) if data['name'] in genes]  # Extract valid genes
    kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Ensure reproducibility
    folds = list(kf.split(valid_genes))  # Get indices for splits
    all_results = []
    k = 100
    for df in np.round(np.arange(0.1,1, 0.1), 1):
        fold_evaluations = useful_functions.run_page_rank_omics(df, graph, folds, valid_genes , k)  
    
        # Collect results for each fold
        for eval_result in fold_evaluations:
            # Add damping factor and fold to the result
            eval_result['damping_factor'] = df
            # Append to the results list
            all_results.append(eval_result)


    results_df = pd.DataFrame(all_results)
    average_results_per_damping_factor = results_df.groupby("damping_factor").mean()
    average_results_per_damping_factor = average_results_per_damping_factor.reset_index()
    average_results_per_damping_factor['n_nodes'] = n_nodes
    average_results_per_damping_factor['n_edges'] = n_edges
    average_results_per_damping_factor['n_tr_genes'] = len(valid_genes)
    average_results_per_damping_factor['correlation_threshold'] = correlation_threshold
    average_results_per_damping_factor.to_csv(out_file, index=False)


def create_out_dir():
    today_date = datetime.today().strftime('%Y_%m_%d')  # Format: YYYY-MM-DD
    omics_output = f"/data/Julia/Omics_output_{today_date}/"
    # Create the folder if it doesn't exist
    if not os.path.exists(omics_output):
        os.makedirs(omics_output)
        print(f"Folder '{omics_output}' created successfully!")
    else:
        print(f"Folder '{omics_output}' already exists.")
    return omics_output

    
def exctract_largest_cc(graph, genes):
    components = list(nx.connected_components(graph))
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
    
    if len(overlap_genes) > 0 :
        print(f"Genes in non-connected components that are also in the training set: {len(overlap_genes)}")
        print(f"These genes are: {overlap_genes}")
    else:
        print("There are no genes in overlap")
        
    largest_cc = max(components, key=len)
    largest_cc_subgraph = graph.subgraph(largest_cc)

    # Verify the size of the largest connected component
    print(f"Largest connected component has {len(largest_cc_subgraph.nodes)} nodes and {len(largest_cc_subgraph.edges)} edges.")
    return largest_cc_subgraph

    
    
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
    omics_output = create_out_dir()
    files = os.listdir(omics_graph_dir)
    for file in files:
        print(f"analysing the graph {file}")
        graph_file = os.path.join(omics_graph_dir, file)
        G = nx.read_graphml(graph_file)
        # Check if the graph is connected (for undirected graphs)
        out_file = f"{omics_output}{os.path.splitext(file)[0]}_angio.csv"
        analyse(G , angiogenic_genes , out_file )
    
        out_file = f"{omics_output}{os.path.splitext(file)[0]}_tip_cell.csv"
        analyse(G , tip_cell_genes , out_file )

if __name__ == "__main__":
    main()
   