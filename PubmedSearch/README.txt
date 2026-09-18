The Pubmed_search_Leo_Bettoni.py script allows you to look for publications in the PubMed database. It is based on a prompt system as the script asks you for several information, including :

    - Would you like to resume the previous search? (In case the connexion to the PubMed servers was lost)
    --> This is useful in case the connexion to the PubMed servers is lost during the search, which can unfortunately happen during long search

    - Do you want to look for GENERAL papers?
    --> If yes, the algorithm will look for publications simply mentioning the query gene's name in its title/abstract

    - How many SPECIFIC topics are you considering? 
    --> You can enter as many keywords as you want, as a result in addition to the general papers the algorithm will look for papers mentioning the query gene's name and your specific keyword (ex: cancer, angiogenesis, diabetes, ...)

    - Would you like to consider ALIASES for the search?
    --> If yes, the algorithm will also consider for each query gene its known aliases (taken from the GeneCards database) to ensure making a thorough search, which however increases the computational runtime
    --> If yes and your input file doesn't list the aliases of the query genes, they will automatically be matched using the Aliases.csv reference file

    - 