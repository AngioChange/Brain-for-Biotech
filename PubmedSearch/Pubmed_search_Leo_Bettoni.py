# 1. Importing libraries
import sys
sys.path.append("../scripts")
import ast
import Pubmed_search_functions as psf
import glob
import os
import pandas as pd
import pyfiglet
import shutil
from Bio import Entrez
from prompt_toolkit import prompt
from time import sleep
from tqdm import tqdm

# 2. Printing a banner
print("\n")
banner = pyfiglet.figlet_format("Pubmed Search")
print(banner, "\n", "(c) Léo Bettoni 2023", "\n"*3)

# 3. Print message error if no DEG list is provided:
if len(sys.argv) != 2:
    sys.exit("Error : missing list of genes file !")
if not os.path.exists(sys.argv[1]):
    sys.exit("Error : file {:s} doesn't exist !".format(sys.argv[1]))

 # 4. Reading the list of DEG and aliases :
input_file = sys.argv[1]
df_pubmed_search = pd.read_csv(input_file)
reference = pd.read_csv("Aliases.csv")

# 5. Getting the name of the file
basename = os.path.basename(input_file)
result_name = os.path.splitext(basename)

# 6. Creating the list of positive answers
answers_yes = ["yes", "YES", "Yes", "y", "Y"]
answers_no = ["no", "NO", "No", "n", "N"]
accepted_answers = answers_yes + answers_no

# 7. Creating the condition list
conditions = []

# 8. Running a previous search that failed ?
sleep(1)
taille = len("Would you like to run/resume the previous search ?")
print(" "*15, "="*(taille+4))
print(" "*15, "* Would you like to run/resume the previous search ?")
print(" "*15, "="*(taille+4), "\n")

retry = input("Your answer : ")
print("\n")

if retry not in accepted_answers:
    while retry not in accepted_answers:
        print("Sorry, I couldn't understand your answer\n")
        retry = input("Your answer : \n")
        os.system("clear")
        print(" "*15, "="*(taille+4), "\n", " "*15,"* Would you like to run/resume the previous search again ?", "\n", " "*15, "="*(taille+4), "\n")


if retry in answers_yes:
    if not os.path.exists("search_parameters.csv"):
        sys.exit("Error : This is your first search, so no previous parameters were saved !")

    else:

        search_parameters = pd.read_csv("search_parameters.csv")
        for parameter in search_parameters.itertuples():

            need_splitting = parameter[1]
            answer_general = parameter[2]
            answer1 = parameter[3]

            if parameter[4] == "none":
                conditions = []
            else:
                conditions = ast.literal_eval(parameter[4])

            answer2 = parameter[5]

            answer3 = parameter[6]
            if answer3 not in answers_yes:
                df_pubmed_search = psf.adding_relevant_aliases(df_pubmed_search, reference)


            threshold = parameter[7]
            email = parameter[8]

        files = glob.glob("parts/part_*.csv")

    if not os.path.exists(f"Pubmed_search_{result_name[0]}"):
        os.mkdir(f"Pubmed_search_{result_name[0]}")


else:

    # 9. Creating a folder for the results if it doesn't exist
    if os.path.exists("Pubmed_search_{}".format(result_name[0])):
        shutil.rmtree("Pubmed_search_{}".format(result_name[0]))
        os.mkdir("Pubmed_search_{}".format(result_name[0]))

    else:
        os.mkdir("Pubmed_search_{}".format(result_name[0]))

    if os.path.exists("parts"):
        shutil.rmtree("parts")


# 10. If there are several input genes the datafile will be splitted
    if len(df_pubmed_search) > 1:
        need_splitting = "yes"
    else:
        need_splitting = "no"

# 11. Storing the different conditions in the list
    sleep(1)

# 12. Need to look for general papers ?
    taille = len("Do you want to look for GENERAL papers ?")
    print(" "*15, "="*(taille+4))
    print(" "*15, "* Do you want to look for GENERAL papers ?")
    print(" "*15, "="*(taille+4), "\n")

    answer_general = input("Your answer : ")
    if answer_general not in accepted_answers:
        while answer_general not in accepted_answers:
            print("Sorry, I couldn't understand your answer\n")
            answer_general = input("Your answer : \n")
            os.system("clear")
            print(" "*15, "="*(taille+4), "\n", " "*15,"* Do you want to look for GENERAL papers ?", "\n", " "*15, "="*(taille+4), "\n")

    sleep(1)

# 13. Need to look for specific condition(s) ?
    print("\n"*3)
    taille = len("How many SPECIFIC topics are you considering ?")
    print(" "*15, "="*(taille+4))
    print(" "*15, "* How many SPECIFIC topics are you considering ? *")
    print(" "*15, "="*(taille+4), "\n")

    answer1 = int(input("Your answer : "))

    for i in range(0, answer1):
        condition = prompt("\n- Topic {} ?\n".format(i+1))
        conditions.append(condition)

    sys.stdout.flush()

# 14. We need to check if aliases should be considered for the search
    sleep(1)

    print("\n"*3)
    taille = len("Would you like to consider ALIASES for the search ?")
    print(" "*15, "="*(taille+4))
    print(" "*15, "* Would you like to consider ALIASES for the search ? *")
    print(" "*15, "="*(taille+4), "\n")

    answer2 = input("Your answer : ")

    if answer2 not in accepted_answers:
        while answer2 not in accepted_answers:
            print("Sorry, I couldn't understand your answer\n")
            answer2 = input("Your answer : \n")
            os.system("clear")
            print(" "*15, "="*(taille+4), "\n", " "*15, "* Would you like to consider ALIASES for the search ?", "\n", " "*15, "="*(taille+4), "\n")

# 15. Adding relevant aliases if necessary
    if answer2 in answers_yes:
        sleep(1)
        print("\n")
        print("- Does your datafile ALREADY have aliases ? \n")
        answer3 = input("Your answer : ")

        if answer3 not in accepted_answers:
            while answer3 not in accepted_answers:
                print("Sorry, I couldn't understand your answer\n")
                answer3 = input("Your answer : \n")
                os.system("clear")
                print(" "*15, "="*(taille+4), "\n", " "*15, "* Does you datafile ALREADY have aliases ?", "\n", " "*15, "="*(taille+4), "\n")

        if answer3 not in answers_yes:
            df_pubmed_search = psf.adding_relevant_aliases(df_pubmed_search, reference)

# 16. Subsplitting the input file if necessary
    if need_splitting in answers_yes:
        files = psf.split_dataframe(df_pubmed_search)


# 17. Selecting a threshold for finding mystery genes
    if answer_general in answers_yes:
        sleep(1)

        print("\n"*3)
        taille = len("How many GENERAL papers to consider for identifying mystery genes ?")
        print(" "*15, "="*(taille+4))
        print(" "*15, "* How many GENERAL papers to consider for identifying mystery genes ? *")
        print(" "*15, "="*(taille+4), "\n")

        threshold = int(input("Your answer : "))


# 18. We need to retrieve the email of the user
    sleep(1)

    print("\n"*3)
    taille = len("What's your email address ? (Mandatory for Pubmed)")
    print(" "*15, "="*(taille+4))
    print(" "*15, "* What's your email address ? (Mandatory for Pubmed) *")
    print(" "*15, "="*(taille+4), "\n")

    email = input("Your answer : ")

    sleep(1)
    print("\n")

# 19. Saving the parameters in a .csv file
    search_parameters = pd.DataFrame()

    search_parameters["Need splitting input file ?"] = [need_splitting]
    search_parameters["Looking for general papers ?"] = [answer_general]
    search_parameters["How many specific topics ?"] = [answer1]
    if answer1 == 0:
        search_parameters["Which ones ?"] = ["none"]
    else:
        search_parameters["Which ones ?"] = [conditions]
    search_parameters["Looking for aliases ?"] = [answer2]

    
    if answer2 in answers_yes:
        search_parameters["Already have aliases ?"] = [answer3]
    else:
        search_parameters["Already have aliases ?"] = ["no"]
    

    if answer_general not in answers_yes:
        search_parameters["Threshold mystery genes ?"] = ["none"]
    else:
        search_parameters["Threshold mystery genes ?"] = [threshold]

    
    search_parameters["email ?"] = [email]

    search_parameters.to_csv("search_parameters.csv", sep=",", index=False)


# 20. Running a Pubmed search for large files
if need_splitting in answers_yes:
    if not os.path.exists("Pubmed_search_{}/subresults".format(result_name[0])):
        os.mkdir("Pubmed_search_{}/subresults".format(result_name[0]))

    for file in tqdm(files, position=0, desc="Processing files"):
        if os.path.exists("Pubmed_search_{}/subresults/part_{}.csv".format(result_name[0], files.index(file) + 1)):
            continue

        else:
            file2 = pd.read_csv(file)

            # 21'. GENERAL PAPERS NEEDED
            if answer_general in answers_yes:
                results = psf.pubmed_search_general(file2, need_splitting, answer2, email)

                # If no specific condition(s) then we can move on to the next file
                if len(conditions) == 0:
                    results.to_csv("Pubmed_search_{}/subresults/part_{}.csv".format(result_name[0], files.index(file) + 1),
                                   sep=",", index=False)

            # If there are specific condition(s) to look for we can carry on with the search
                else:
                    results2 = psf.pubmed_search_conditions(file2, need_splitting, answer2, conditions, email)
                    results3 = pd.merge(results, results2, on="genes")
                    results3.to_csv("Pubmed_search_{}/subresults/part_{}.csv".format(result_name[0], files.index(file) + 1),
                                    sep=",", index=False)

                # 21''. GENERAL PAPERS NOT NEEDED
            else:
                if len(conditions) == 0:
                    sys.exit("==> No general papers needed and no specific topic(s) provided !")

                else:
                    results = psf.pubmed_search_conditions(file2, need_splitting, answer2, conditions, email)
                    results.to_csv("Pubmed_search_{}/subresults/part_{}.csv".format(result_name[0], files.index(file) + 1),
                                   sep=",", index=False)

    # Merging the subresults and removing all intermediary files/folders
    result_files = glob.glob("Pubmed_search_{}/subresults/*.csv".format(result_name[0]))
    df = pd.concat(map(pd.read_csv, result_files), ignore_index=True)
    df.to_csv("Pubmed_search_{}/results_pubmed_search_{}.csv".format(result_name[0], result_name[0]),
              sep=",", index=False)

    # Applying a threshold to select potential mystery genes
    mystery_genes = df[df["Publications general"] <= threshold]
    mystery_genes.to_csv("Pubmed_search_{}/results_pubmed_search_{}_mystery_genes_threshold_{}_papers.csv"
                         .format(result_name[0], result_name[0], threshold), index=False, sep=",")

    # Removing all intermediary files and folders
    shutil.rmtree("parts")
    shutil.rmtree("Pubmed_search_{}/subresults".format(result_name[0]))

    # Exiting the program
    sys.exit("==> PUBMED SEARCH COMPLETE !")

# 22. Running a Pubmed search for short files
elif need_splitting not in answers_yes:

    # 22'. GENERAL PAPERS NEEDED
    if answer_general in answers_yes:
        results = psf.pubmed_search_general(df_pubmed_search, need_splitting, answer2, email)

        # If no specific condition(s) then we can save this file and end the program
        if len(conditions) == 0:
            results.to_csv("Pubmed_search_{}/results_pubmed_search_{}.csv".format(result_name[0], result_name[0]), sep=",",
                           index=False)

            # Applying a threshold to select potential mystery genes
            mystery_genes = results[results["Publications general"] <= threshold]
            mystery_genes.to_csv("Pubmed_search_{}/results_pubmed_search_{}_mystery_genes_threshold_{}_papers.csv".format(result_name[0], result_name[0], threshold), index=False, sep=",")

            # Exiting the program
            sys.exit("==> PUBMED SEARCH COMPLETE !")

        else:
            results2 = psf.pubmed_search_conditions(df_pubmed_search, need_splitting, answer2, conditions, email)
            results3 = pd.merge(results, results2, on="genes")
            results3.to_csv("Pubmed_search_{}/results_pubmed_search_{}.csv".format(result_name[0], result_name[0]),
                            sep=",", index=False)

            # Applying a threshold to select potential mystery genes
            mystery_genes = results3[results3["Publications general"] <= threshold]
            mystery_genes.to_csv("Pubmed_search_{}/results_pubmed_search_{}_mystery_genes_threshold_{}_papers.csv"
                                 .format(result_name[0], result_name[0], threshold), sep=",", index=False)

            # Exiting the program
            sys.exit("==> PUBMED SEARCH COMPLETE !")

    else:
        # 22''. GENERAL PAPERS NOT NEEDED
        if len(conditions) == 0:
            sys.exit("==> No general papers needed and no specific topic(s) provided !")

        else:
            results = psf.pubmed_search_conditions(df_pubmed_search, need_splitting, answer2, conditions, email)
            results.to_csv("Pubmed_search_{}/results_pubmed_search_{}.csv".format(result_name[0], result_name[0]), sep=",", index=False)

            # Exiting the program
            sys.exit("==> PUBMED SEARCH COMPLETE !")
