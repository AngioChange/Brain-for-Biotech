"""Module to run a nice Pubmed search"""

import ast
import glob
import os
import sys
import pandas as pd
import numpy as np

from prompt_toolkit import prompt
from Bio import Entrez
from tqdm import tqdm
from sympy import symbols, Eq, solve

DATE = 17112023


def adding_relevant_aliases(dataframe, aliases_file):
    """Adds relevant aliases to the input genes if necessary"""
    if dataframe.columns[0] != "genes": ##to make sure the column name is the same as the aliases df
            dataframe = dataframe.rename(columns = {dataframe.columns[0] : "genes"}) 
        
    dataframe[dataframe.columns[0]] = dataframe[dataframe.columns[0]].str.upper() ##to merge with aliases genes need to be in uppercase
    dataframe_merged = dataframe.merge(aliases_file, on='genes', how='inner')
    dataframe_merged = dataframe_merged[["genes", "aliases"]]
    dataframe_merged = dataframe_merged.drop_duplicates()

    #if some genes were not merged in the process (because they are absent from the aliase.csv file), we'll manually add them without aliases
    if len(dataframe_merged.columns[0]) != len(dataframe.columns[0]) + len(aliases_file.columns[0]):
        aliases = dataframe_merged.genes.unique()
        missing_aliases_IDs = dataframe[~dataframe["genes"].isin(aliases)]

        for gene in missing_aliases_IDs.itertuples():
            df_add = pd.DataFrame({"genes": gene[1], "aliases": "NO ALIAS"}, index = [0])
            dataframe_merged = pd.concat([dataframe_merged, df_add], ignore_index = True)
    
    return dataframe_merged


def pubmed_search_general(dataframe, need_splitting, answer_alias, email):
    """Runs a Pubmed search and gives GENERAL papers""" 
    answers_yes = ["yes", "YES", "Yes", "y", "Y"]
    genes = []
    aliases = []
    IDs_list = []

    if need_splitting in answers_yes:
        pos_arg = 1
    else:
        pos_arg = 0

    for gene in tqdm(dataframe.itertuples(), position = pos_arg, desc = "Processing the input gene(s) for GENERAL publications", leave = False):
        genes.append(gene[1])
        
        if answer_alias in answers_yes:
            if gene[2] == "NO ALIAS":
                aliases.append(gene[2])
            else:
                aliases.append(ast.literal_eval(gene[2]))
            
        IDs = []
        Entrez.email = email
        handle_search = Entrez.esearch(db="pubmed", sort="Relevance", term='"{}"[Title/Abstract]'.format(gene[1]), usehistory="y", retmax = 10000)
        record_search = Entrez.read(handle_search)

        for record in record_search['IdList']:
            IDs.append(record)

    
        # Checking for aliases
        if answer_alias not in answers_yes:
            IDs_list.append(IDs)
    
        else:

            if gene[2] == "NO ALIAS":
                IDs_list.append(IDs)
        
            else:

                for alias in tqdm(ast.literal_eval(gene[2]), position = pos_arg + 1, desc = "Processing aliases", leave = False):
                    Entrez.email = email
                    handle_search = Entrez.esearch(db="pubmed", sort="Relevance", term='"{}"[Title/Abstract]'.format(alias), usehistory="y", retmax = 10000)
                    record_search = Entrez.read(handle_search)

                    for record in record_search['IdList']:
                        if record not in IDs:
                            IDs.append(record)

                IDs_list.append(IDs)

    #Removing potential duplicates in the results and computing the actual number of unique publications
    number_of_publications = []
    IDs_list_unique = []

    for sublist in IDs_list:
        sublist = list(dict.fromkeys(sublist))
        IDs_list_unique.append(sublist)
        number_of_publications.append(len(sublist))

    results = pd.DataFrame()
    results["genes"] = genes

    if answer_alias in answers_yes:
        results["aliases"] = aliases
        results["Publications general"] = number_of_publications
        results["Pubmed IDs general"] = IDs_list_unique

        return results

    else:
        results["Publications general"] = number_of_publications
        results["Pubmed IDs general"] = IDs_list_unique

        return results


def pubmed_search_conditions(dataframe, need_splitting, answer_alias, conditions, email):
    """Runs a pubmed search for specific condition(s)"""
    answers_yes = ["yes", "YES", "Yes", "y", "Y"]
    genes = []
    aliases = []

    if need_splitting in answers_yes:
        pos_arg = 1
    else:
        pos_arg = 0

    for gene in dataframe.itertuples():
        genes.append(gene[1])

    results = pd.DataFrame()
    results["genes"] = genes
    
    for condition in conditions:
        IDs_list = []

        for gene in tqdm(dataframe.itertuples(), position = pos_arg, desc = "Processing the input gene(s) for {} publications".format(condition.upper()), leave = False):
            Entrez.email = email
            IDs = []
            handle_search = Entrez.esearch(db="pubmed", sort="Relevance", term='"{}"[Title/Abstract]+"{}"[Title/Abstract]'.format(gene[1], condition), \
                usehistory="y", retmax = 10000)
            record_search = Entrez.read(handle_search)
            record_search["IdList"] = list(dict.fromkeys(record_search["IdList"]))

            for record in record_search['IdList']:
                IDs.append(record)

            

            # Checking aliases
            if answer_alias not in answers_yes:
                IDs_list.append(IDs)
            
            else:

                if gene[2] == "NO ALIAS":
                    IDs_list.append(IDs)
        
                else:
                    for alias in tqdm(list(ast.literal_eval(gene[2])), position = pos_arg + 1, desc = "Processing aliases", leave = False):
                        Entrez.email = email
                        handle_search = Entrez.esearch(db="pubmed", sort="Relevance", term='"{}"[Title/Abstract]'.format(alias), usehistory="y", retmax = 10000)
                        record_search = Entrez.read(handle_search)
                        record_search["IdList"] = list(dict.fromkeys(record_search["IdList"]))


                        for record in record_search['IdList']:
                            IDs.append(record)
    

                    IDs_list.append(IDs)

        #Removing potential duplicates in the sublist of PMIDs
        number_of_publications = []
        IDs_list_unique = []

        for sublist in IDs_list:
            sublist = list(dict.fromkeys(sublist))
            IDs_list_unique.append(sublist)
            number_of_publications.append(len(sublist))
        
        results["Publications in {}".format(condition)] = number_of_publications
        results["Pubmed IDs {}".format(condition)] = IDs_list_unique
        
    return results


    
def split_dataframe(dataframe):
    """Splits a dataframe into sub-files"""
    #Creating the folder for subfiles
    directory = os.getcwd()
    os.mkdir("{}/parts".format(directory))
    
    #Splitting the file into subfiles
    df_split = np.array_split(dataframe, len(dataframe))

    #Saving the results
    for i in range(0, len(dataframe)):
        df = df_split[i]
        df = pd.DataFrame(df, columns = ["genes", "aliases"])
        df.to_csv("{}/parts/part_{}.csv".format(directory, (i + 1)), sep = ",", index = False)
        
    files = glob.glob("{}/parts/*.csv".format(directory))
    
    return files





    

