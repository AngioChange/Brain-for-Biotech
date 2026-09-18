import argparse
import json
import gc
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import requests
import seaborn as sns
import time
from matplotlib.ticker import MaxNLocator
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
from tqdm import tqdm



def evaluate_fold(test_genes, ppr_scores_with_gene_names, k=100):
    # 1. Sort genes by PPR scores
    sorted_genes = sorted(ppr_scores_with_gene_names.items(), key=lambda x: x[1], reverse=True)
    ranked_genes = [gene for gene, score in sorted_genes]
    
    # 2. Compute ranks of test genes
    ranks = [ranked_genes.index(gene) + 1 for gene in test_genes if gene in ranked_genes]
    
    # Mean and median rank
    mean_rank = sum(ranks) / len(ranks) if ranks else None
    median_rank = sorted(ranks)[len(ranks) // 2] if ranks else None
    
    # Precision@k and Recall@k
    top_k_genes = ranked_genes[:k]
    precision_at_k = len([gene for gene in test_genes if gene in top_k_genes]) / k
    recall_at_k = len([gene for gene in test_genes if gene in top_k_genes]) / len(test_genes)
    
    # AUC-ROC and AUC-PR
    y_true = [1 if gene in test_genes else 0 for gene in ranked_genes]
    y_scores = [ppr_scores_with_gene_names[gene] for gene in ranked_genes]
    auc_roc = roc_auc_score(y_true, y_scores) if any(y_true) else None
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    auc_pr = auc(recall, precision) if any(y_true) else None
    
    return {
        "mean_rank": mean_rank,
        "median_rank": median_rank,
        "precision_at_k": precision_at_k,
        "recall_at_k": recall_at_k,
        "auc_roc": auc_roc,
        "auc_pr": auc_pr
    }

def run_page_rank(df, graph , k, folds, valid_genes):
    fold_evaluations = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        fold = fold_idx + 1
        print(f"Processing Fold {fold}...")
    
        train_genes = [valid_genes[i] for i in train_idx]  # Training set (seed nodes)
        test_genes = [valid_genes[i] for i in test_idx]    # Test set
    
        # Build the personalization vector (1 for training genes, 0 for others)
        personalization = {node: 0 for node in graph.nodes}
        for node, data in graph.nodes(data=True):
            if 'name' in data and data['name'] in train_genes:
                personalization[node] = 1.0
    
        # Normalize the personalization vector
        total = sum(personalization.values())
        personalization = {k: v / total for k, v in personalization.items()}
    
        # Run Personalized PageRank
        ppr_scores = nx.pagerank(graph, alpha=df, personalization=personalization)
    
        ppr_scores_with_gene_names = {
            graph.nodes[node]['name']: score
            for node, score in ppr_scores.items()
            if 'name' in graph.nodes[node] and graph.nodes[node]['name'] not in train_genes # Ensure the node has a 'name' attribute
        }
        
        evaluation = evaluate_fold(test_genes, ppr_scores_with_gene_names, k=k)
        fold_evaluations.append({"fold": fold, **evaluation})
    return fold_evaluations




def run_page_rank_omics(df, graph , folds, valid_genes , k = 100):
    fold_evaluations = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        fold = fold_idx + 1
        
        train_genes = [valid_genes[i] for i in train_idx]  # Training set (seed nodes)
        test_genes = [valid_genes[i] for i in test_idx]    # Test set
    
        # Build the personalization vector (1 for training genes, 0 for others)
        personalization = {node: 0 for node in graph.nodes}
        for node, data in graph.nodes(data=True):
            if 'name' in data and data['name'] in train_genes:
                personalization[node] = 1.0
    
        # Normalize the personalization vector
        total = sum(personalization.values())
        personalization = {k: v / total for k, v in personalization.items()}
    
        # Run Personalized PageRank
        ppr_scores = nx.pagerank(graph, alpha=df, personalization=personalization)
    
         # remove training genes from the page rank scores
        ppr_scores_filtered = {
            graph.nodes[node]['name']: score
            for node, score in ppr_scores.items()
            if 'name' in graph.nodes[node] and graph.nodes[node]['name'] not in train_genes 
        }
       
        evaluation = evaluate_fold(test_genes, ppr_scores_filtered, k=k)
        fold_evaluations.append({"fold": fold, **evaluation})
    return fold_evaluations

def run_page_rank_databases(df, graph ,folds, valid_genes, k=100):
    fold_evaluations = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        fold = fold_idx + 1
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
            gene: score
            for gene, score in ppr_scores.items()
            if gene not in train_genes # Ensure the node not in the training genes
        }
        evaluation = evaluate_fold(test_genes, ppr_scores_filtered, k=k)
        fold_evaluations.append({"fold": fold, **evaluation})
    return fold_evaluations



def predict_all(df , graph , valid_genes):
    personalization = {node: 0 for node in graph.nodes}
    for node in graph.nodes(data=False):
        if node in valid_genes:
            personalization[node] = 1.0
    
    # Normalize the personalization vector
    total = sum(personalization.values())
    personalization = {k: v / total for k, v in personalization.items()}
    
    # Run Personalized PageRank
    ppr_scores = nx.pagerank(graph, alpha=df, personalization=personalization , weight="weight")

    df = pd.DataFrame(list(ppr_scores.items()), columns=["Gene", "Score"])
    df = df.sort_values(by="Score", ascending=False)
    
    df["Label"] = df["Gene"].apply(lambda gene: 1 if gene in valid_genes else 0)
    
    return df



def merge(ppi_graph, omics_graph):
    merged_graph = nx.Graph()  # Use DiGraph() if your graph is directed

    # Add edges and weights from the BioGRID graph
    for u, v, data in ppi_graph.edges(data=True):
        weight = data.get('weight')  # Default weight is 1 if not present
        if merged_graph.has_edge(u, v):
            merged_graph[u][v]['weight'] += weight  # Combine weights
        else:
            merged_graph.add_edge(u, v, weight=weight)
    
    for u, v  in omics_graph.edges(data=False):
        weight = 1  # Default weight is 1 if not present
        u = omics_graph.nodes[u]['name']
        v = omics_graph.nodes[v]['name']
        if merged_graph.has_edge(u, v):
            merged_graph[u][v]['weight'] += weight  # Combine weights
        else:
            merged_graph.add_edge(u, v, weight=weight)
            
    return merged_graph


def permutations_predict_all(interest, out_dir, merged_graph):
    genes = list(merged_graph.nodes())
    shuffled_genes = genes[:]
    print("Predicting in real graph real genes")

    predictions = predict(alpha, merged_graph, genes_of_interest[interest])
    out_file = f"{out_dir}/PageRank_predictions_real_graph_real_genes_{interest}.csv"
    predictions.to_csv(out_file, index = False)

    for i in range(10):
        print(f"Running permutation {i+1}")
        random.shuffle(shuffled_genes)

        # Relabel nodes with shuffled gene names
        shuffled_graph = merged_graph.copy()

        # Try edge swapping with decreasing swap numbers until successful
        max_swaps = len(merged_graph.edges())
        for swap_factor in [5, 3, 2, 1]:
            try:
                nswap = max_swap * swap_factor
                nx.double_edge_swap(shuffle_graph, nswap = nswap, max_tries = nswap*10)
                print(f"Successfully swapped with factor {swap_factor}")
                break
            except nx.NetworkXError:
                print(f"Swap factor {swap_factor} failed, trying smaller factor ...")
                continue
            else:
                print("Warning:  Could not perform edge swapping, using original graph structure")

            sampled_genes = random.sample(shuffled_genes, n_sample_genes[interest])

            print("Predicting in shuffled graph shuffled genes")
            predictions = predict(alpha, shuffled_graph, sampled_genes)
            out_file = f"{out_dir}/PageRank_predictions_perm_graph_perm_genes_{i+1}_{interest}.csv"
            predictions.to_csv(out_file, index = False)

            print("Predicting in real graph shuffled genes")
            predictions = predict(alpha, merged_graph, sampled_genes)
            out_file = f"{out_dir}/PageRank_predictions_real_graph_perm_genes_{i+1}_{interest}.csv"
            predictions.to_csv(out_file, index = False)

            print("Predicting in shuffled graph real genes")
            predictions = predict(alpha, shuffled_graph, genes_of_interest[interest])
            out_file = f"{out_dir}/PageRank_predictions_perm_graph_real_genes_{i+1}_{interest}.csv"
            predictions.to_csv(out_file, index = False)

            del shuffled_graph
            gc.collect()
    



def permutations_predict_test(process, alpha, out_dir, merged_graph, valid_genes):
    genes = list(merged_graph.nodes())
    shuffled_genes = genes[:]

    for i in range(10):
        print(f"Running permutation {i+1}")

        predictions = predict_in_test(alpha, merged_graph, valid_genes)
        out_file = f"{out_dir}/Real_data/PageRank_predictions_in_test_real_graph_real_genes_{i+1}_{process}.csv"
        predictions.to_csv(out_file, index = False)

        # Relabel nodes with shuffled gene names
        shuffled_graph = merged_graph.copy()

        # Try edge swapping with decreasing swap numbers until it's successful
        max_swaps = len(merged_graph.edges())
        for swap_factor in [5, 3, 2, 1]:
            try:
                nswap = max_swaps * swap_factor
                nx.double_edge_swap(shuffled_graph, nswap = nswap, max_tries = nswap * 10)
                print(f"Successfully swapped with factor {swap_factor}")
                break
            except nx.NetworkXError:
                print(f"Swap factor {swap_factor} failed, trying smaller factor ...")
                continue
        else:
            print("Warning: Could not perform edge swapping, using original graph structure")

        sampled_genes = random.sample(shuffled_genes, len(valid_genes))

        print("Predicting in shuffled graph shuffled genes")
        predictions = predict_in_test(alpha, shuffled_graph, sampled_genes)
        out_file = f"{out_dir}/Perm_graph_perm_genes/PageRank_predictions_in_test_perm_graph_perm_genes_{i+1}_{process}.csv"
        predictions.to_csv(out_file, index = False)

        print("Predicting in real graph shuffled genes")
        predictions = predict_in_test(alpha, merged_graph, sampled_genes)
        out_file = f"{out_dir}/Real_graph_perm_genes/PageRank_predictions_in_test_real_graph_perm_genes_{i+1}_{process}.csv"
        predictions.to_csv(out_file, index = False)

        print("Predicting in shuffled graph real genes")
        predictions = predict_in_test(alpha, shuffled_graph, valid_genes)
        out_file = f"{out_dir}/Perm_graph_real_genes/PageRank_predictions_in_test_perm_graph_real_genes_{i+1}_{process}.csv"
        predictions.to_csv(out_file, index = False)

        del shuffled_graph
        gc.collect()


    
def predict_in_test(df, graph, valid_genes):
    # split valid genes into 80/20 train/test split
    # first copy valid genes and permute them
    shuffled_genes = valid_genes[:]
    np.random.shuffle(shuffled_genes)
    train_genes = shuffled_genes[:int(len(shuffled_genes) * 0.8)]
    test_genes = shuffled_genes[int(len(shuffled_genes) * 0.8):]
    
    # run PageRank on the graph with the train genes as personalization
    personalization = {node: 0 for node in graph.nodes}
    for node in graph.nodes(data = False):
        if node in train_genes:
            personalization[node] = 1.0

    # Normalize the personalization vector
    total = sum(personalization.values())
    personalization = {k: v / total for k, v in personalization.items()}

    # Run the PPR algorithm
    ppr_scores = nx.pagerank(graph, alpha = df, personalization = personalization, weight = "weight")
                             
    # Filter the scores to only include test genes
    ppr_scores_filtered = {
        gene: score
        for gene, score in ppr_scores.items()
        if gene in test_genes # Ensure the node is in the test genes
    }

    # Create a DataFrame with the scores
    df = pd.DataFrame(list(ppr_scores_filtered.items()), columns = ["Gene", "Score"])
    df = df.sort_values(by = "Score", ascending = False)
    df["Label"] = df["Gene"].apply(lambda gene: 1 if gene in test_genes else 0)
    return df
            


def analyse(merged_graph, genes):
    n_nodes = len(merged_graph.nodes())
    n_edges = len(merged_graph.edges())
    
    valid_genes = [node for node in merged_graph.nodes(data=False) if node in genes] 
    kf = KFold(n_splits=5, shuffle=True, random_state=42)  # Ensure reproducibility
    folds = list(kf.split(valid_genes))  # Get indices for splits

    all_results = []
    k = 100
    damping_factors = np.round(np.arange(0.1, 1, 0.1), 1) 
    for df in tqdm(damping_factors, desc = "Applying damping factors ..."):
        fold_evaluations = run_page_rank_databases(df, merged_graph, folds, valid_genes , k)  # Function 
    
        # Collect results for each fold
        for eval_result in fold_evaluations:
            # Add damping factor and fold to the result
            eval_result['damping_factor'] = df
            # Append to the results list
            all_results.append(eval_result)

    results_df = pd.DataFrame(all_results)
    average_results_per_damping_factor = results_df.groupby("damping_factor").mean()
    average_results_per_damping_factor['n_nodes'] = n_nodes
    average_results_per_damping_factor['n_edges'] = n_edges
    average_results_per_damping_factor['n_tr_genes'] = len(valid_genes)
    average_results_per_damping_factor = average_results_per_damping_factor.sort_values(by=[ "auc_roc"], ascending=[False])
    # Reset index for cleaner display
    average_results_per_damping_factor = average_results_per_damping_factor.reset_index()
    return average_results_per_damping_factor



# Enrichr API Function for Manhattan Plot and Bar Chart
# Takes a gene list and Enrichr libraries as input
# Finds ALL significantly enriched pathways
def Enrichr_API(enrichr_gene_list, all_libraries):

    all_terms = []
    all_pvalues =[] 
    all_adjusted_pvalues = []

    for library_name in all_libraries : 
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(enrichr_gene_list)
        description = ''
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)
        time.sleep(0.5)
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        user_list_id = data['userListId']
        short_id = data["shortId"]
        gene_set_library = library_name
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)

        if len(data[library_name]) == 0:
            raise NoResults

        short_results_df  = pd.DataFrame(data[library_name])
        all_terms_list = []
        all_pvalues_list = []
        all_adjusted_pvalues_list = []

        all_terms = []
        all_pvalues = [] 
        all_adjusted_pvalues = []

        for value in short_results_df.itertuples():
            if value[3] < 0.05:
                all_terms_list.append(value[2])
                all_pvalues_list.append(value[3])
                all_adjusted_pvalues_list.append(value[7])

        all_terms.append(all_terms_list)
        all_pvalues.append(all_pvalues_list)
        all_adjusted_pvalues.append(all_adjusted_pvalues_list)
        
        results_df  = pd.DataFrame(data[library_name])
        # adds library name to the data frame so the libraries can be distinguished
        results_df['library'] = library_name.replace('_', '')

    return [results_df, short_results_df, all_terms, all_pvalues, all_adjusted_pvalues, str(short_id)]




# Enrichr API Function for Manhattan Plot and Bar Chart
# Takes a gene list and Enrichr libraries as input
# Finds top n significantlty enriched pathways
def Enrichr_API_top_n(enrichr_gene_list, all_libraries, threshold_pathways):

    all_terms = []
    all_pvalues =[] 
    all_adjusted_pvalues = []

    for library_name in all_libraries : 
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(enrichr_gene_list)
        description = ''
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)
        time.sleep(0.5)
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        user_list_id = data['userListId']
        short_id = data["shortId"]
        gene_set_library = library_name
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise APIFailure

        data = json.loads(response.text)

        if len(data[library_name]) == 0:
            raise NoResults

        short_results_df  = pd.DataFrame(data[library_name][0:threshold_pathways])
        all_terms_list = []
        all_pvalues_list = []
        all_adjusted_pvalues_list = []

        all_terms = []
        all_pvalues = [] 
        all_adjusted_pvalues = []

        for value in short_results_df.itertuples():
            if value[3] < 0.05:
                all_terms_list.append(value[2])
                all_pvalues_list.append(value[3])
                all_adjusted_pvalues_list.append(value[7])

        all_terms.append(all_terms_list)
        all_pvalues.append(all_pvalues_list)
        all_adjusted_pvalues.append(all_adjusted_pvalues_list)
        
        results_df  = pd.DataFrame(data[library_name])
        # adds library name to the data frame so the libraries can be distinguished
        results_df['library'] = library_name.replace('_', '')

    return [results_df, short_results_df, all_terms, all_pvalues, all_adjusted_pvalues, str(short_id)]



# Bar Chart Functions
# Takes all terms, all p-values, all adjusted p-values, plot title, Enrichr libraries, and specified figure format
# Displays ALL significantly enriched pathways
def enrichr_figure(all_terms, all_pvalues, all_adjusted_pvalues, plot_name, all_libraries, bar_color, scope, x_param, y_param): 
    # Bar colors
    if bar_color != 'lightgrey':
        bar_color_not_sig = 'lightgrey'
        edgecolor=None
        linewidth=0
        bar_color_not_sig = 'white'
        edgecolor='black'
        linewidth=1    

    plt.figure(figsize=(x_param, y_param))
    
    i = 0
    bar_colors = [bar_color if (x < 0.05) else bar_color_not_sig for x in all_pvalues[i]]
    fig = sns.barplot(x=np.log10(all_pvalues[i])*-1, y=all_terms[i], palette=bar_colors, edgecolor=edgecolor, linewidth=linewidth)
    fig.axes.get_yaxis().set_visible(False)
    fig.set_title(f"Top {scope} significantly enriched pathways - {all_libraries[i].replace('_', ' ')}", fontsize=26)
    fig.set_xlabel('−log₁₀(p‐value)', fontsize=25)
    fig.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tick_params(axis='x', which='major', labelsize=20)
    if max(np.log10(all_pvalues[i])*-1)<1:
        fig.xaxis.set_ticks(np.arange(0, max(np.log10(all_pvalues[i])*-1), 0.1))
    for ii,annot in enumerate(all_terms[i]):
        if all_adjusted_pvalues[i][ii] < 0.05:
            annot = '  *'.join([annot, str(str(np.format_float_scientific(all_pvalues[i][ii], precision=2)))]) 
        else:
            annot = '  '.join([annot, str(str(np.format_float_scientific(all_pvalues[i][ii], precision=2)))])

        title_start= max(fig.axes.get_xlim())/200
        fig.text(title_start, ii, annot, ha='left', wrap = True, fontsize = 26)

    fig.spines['right'].set_visible(False)
    fig.spines['top'].set_visible(False)

    
    plt.savefig(plot_name, bbox_inches = 'tight')
    
    # Show plot 
    plt.show()  


def compute_metrics(predictions, labels):
    predictions_cls = (predictions > 0.5).astype(int)
    accuracy = np.mean(predictions_cls == labels)
    roc_auc = roc_auc_score(labels, predictions)
    precision_recall = precision_recall_curve(labels, predictions)
    pr_auc = auc(precision_recall[1], precision_recall[0])
    return accuracy, roc_auc, pr_auc


def graph_to_matrix(graph):
    genes = list(graph.nodes())
    matrix = np.zeros((len(genes), len(genes)))
    gene_index = {gene: idx for idx, gene in enumerate(genes)}

    for u, v, data in graph.edges(data = True):
        if "weight" in data:
            weight = data["weight"]
        else:
            weight = 1.0 # Default weight if not specified
        matrix[gene_index[u], gene_index[v]] = weight
        matrix[gene_index[v], gene_index[u]] = weight # Assuming undirected graph

    return pd.DataFrame(matrix, index = genes, columns = genes)