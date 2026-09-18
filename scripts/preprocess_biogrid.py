import os
import networkx as nx
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm
import json

annot_graph_dir = "G:/Leo/BfBio_reloaded/data/graphs/annotation-based/BioGrid/"
annot_network_dir = "G:/Leo/BfBio_reloaded/data/annotation-based/BioGrid/"

# Step 1: Remove NaNs and negative scores from numeric methods 
def filter_invalid_and_nan_scores(df, target_systems):
    """
    Removes rows with NaN scores for the specified target systems
    and removes rows where the Score <= 0 for all systems.
    """
    df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
    # Remove rows with NaN scores for the target systems
    df = df[~((df['Experimental System'].isin(target_systems)) & (df['Score'].isna()))]
    print(df.shape)
    # Remove rows where Score <= 0 for all systems
    df = df[~((df['Experimental System'].isin(target_systems)) & (df['Score'] <= 0)) ]
    print(df.shape)
    return df

# Step 2: Normalize scores per experimental system
def normalize_scores(df, target_systems):
    # normalize the score per target system, because the range differs per system
    # bring everething on the same rage
    max_score_overall = df[df['Experimental System'].isin(target_systems)]['Score'].max()
    print(max_score_overall)
    df['Normalized Score'] = df.groupby('Experimental System')['Score'].transform(
        lambda x: (x / x.max() * max_score_overall) if x.max() > 0 else x
    )
    return df

# Step 3: Assign most frequent values to ['Cross-Linking-MS (XL-MS)', 'Co-crystal Structure']
# Since there are not score for these experimental systems
def assign_most_frequent_scores(df, source_systems, target_systems):
    # Find the most frequent values in the source systems
    most_frequent_scores = df[df['Experimental System'].isin(source_systems)]['Normalized Score'].mode()
    if not most_frequent_scores.empty:
        most_frequent_value = most_frequent_scores.iloc[0]  # Use the first mode if multiple modes exist
        # Assign the most frequent value to rows in the target systems
        df.loc[df['Experimental System'].isin(target_systems), 'Normalized Score'] = most_frequent_value
    return df
        
def create_graph():
    biogrid = pd.read_csv(annot_network_dir + "BIOGRID-ALL-4.4.246.tab3.zip" , sep='\t')
    print(f"The original graph has the shape: {biogrid.shape}")
    biogrid = biogrid[(biogrid['Organism ID Interactor A'] == 9606) & (biogrid['Organism ID Interactor B'] == 9606)]
    print(f"After filtering for 9606 and methods graph has the shape: {biogrid.shape}")
    
    trustworthy_methods = [
    'Affinity Capture-MS', 'Co-crystal Structure', 
    'Cross-Linking-MS (XL-MS)', 'Proximity Label-MS',
    'Co-purification','Two-hybrid', 
    'Synthetic Lethality',
    'Proximity Label-MS',
    'FRET'
    ]
    biogrid_filtered = biogrid[
        (biogrid['Experimental System'].isin(trustworthy_methods)) &
        (biogrid['Experimental System Type'].isin(['physical' , 'genetic']) )
    ]
    print(f"After filtering for trustworthy methods, the dataset shape was: " , biogrid_filtered.shape)
    
    methods = set(biogrid_filtered['Experimental System'].to_list())
    
    biogrid_filtered['Score'] = biogrid_filtered['Score'].replace('-', np.nan)
    cols = ['Official Symbol Interactor A','Official Symbol Interactor B' , 'Score', 'Experimental System' ]
    biogrid_filtered = biogrid_filtered[cols]
    
    method_type = {}
    for method in methods:
        scores = biogrid_filtered[biogrid_filtered['Experimental System'] == method]['Score']
        if sum(scores.isna()) == len(scores) :
            method_type[method] = 'not_numeric'
        else:
            method_type[method] = 'numeric'
            
    numeric_methods = [key for key, value in method_type.items() if value == 'numeric']
    not_numeric_methods = [key for key, value in method_type.items() if value == 'not_numeric']
     
    source_systems = numeric_methods
    target_systems = not_numeric_methods

    # Apply the steps
    biogrid_filtered = filter_invalid_and_nan_scores(biogrid_filtered, source_systems)            # Step 1
    biogrid_filtered = normalize_scores(biogrid_filtered, source_systems)             # Step 2
    biogrid_filtered = assign_most_frequent_scores(biogrid_filtered, source_systems, target_systems)  # Step 3
    
    print(biogrid_filtered.head() )
    biogrid_filtered['Normalized Score'] = np.log10(biogrid_filtered['Normalized Score'] + 1 ) 
    'Official Symbol Interactor A','Official Symbol Interactor B' , 'Score', 'Experimental System'
    cols = ['Official Symbol Interactor A','Official Symbol Interactor B' , 'Normalized Score']
    biogrid_filtered = biogrid_filtered[cols]
    biogrid_filtered.rename(columns={"Official Symbol Interactor A": "source", "Official Symbol Interactor B":
                                     "target", "Normalized Score" : "confidence"}, 
                            inplace=True)
    
    biogrid_filtered = biogrid_filtered.sort_values(by="confidence", ascending=False)

    # Step 2: Remove duplicates based on 'gene_name', keeping the first occurrence
    biogrid_filtered = biogrid_filtered.drop_duplicates(subset=["source" , "target"], keep="first")

    # Remove self-interactions
    biogrid_filtered =  biogrid_filtered[biogrid_filtered["source"] != biogrid_filtered["target"]]
    print(f"The shape of the data is {biogrid_filtered.shape}")
    print(f" Min score after normalization is: {biogrid_filtered['confidence'].min()}")
    print(f" Max score after normalization is: {biogrid_filtered['confidence'].max()}")
    
    biogrid_filtered.to_csv(f"{annot_network_dir}biogrid_network.csv" , index = False)
