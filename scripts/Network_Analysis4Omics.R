library(Matrix)
library(patchwork)
library(plyr)
library(Seurat)
library(stringr)
library(tidyselect)
library(igraph)
library(pROC)
library(dplyr)
library(PRROC)
library(Metrics)
library(MLmetrics) 



# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Parse the arguments
if (length(args) < 2) {
  stop("You must provide at least 2 arguments: correlation threshold, interest,  (and graph optional if you want to write graph to file)")
}

corr_threshold <- as.numeric(args[1])  # Convert the first argument to numeric
interest <- args[2]                    # Keep the fourth argument as string (or convert if needed)

write_graph = FALSE
if(args[3] == 'graph') write_graph = TRUE
n_folds = 5

# Print the arguments (optional, for debugging)
cat("Correlation threshold:", corr_threshold, "\n")
cat("Number of Folds:", n_folds, "\n")
cat("Interest:", interest, "\n")
cat("Write graph:", write_graph, "\n")
# Your script logic here
# Example: Run the PageRank algorithm with these arguments
cat("Running PageRank with the provided arguments...\n")


setwd("~/BfBio/")
datasets_dir<- "datasets"
files<-list.files(datasets_dir)
studies<-c()

df_symbols <- read.csv("TrainingGenes/Symbols.csv")

graph_outdir<-"graphs/omics/"
# create the graph folder only in case the graph parameter is TRUE
if (write_graph){
  if (!file.exists(graph_outdir)){
    dir.create(graph_outdir)
    print(paste0(graph_outdir , " is created"))
  }
  
}

outdir<-paste0("graphs/omics/", interest, "/PR_Analysis")

if (!file.exists(outdir)){
  dir.create(outdir)
  print(paste0(outdir , " is created"))
}


stalk_cell_genes <- read.csv("TrainingGenes/training_genes_stalk_cell.csv")
stalk_cell_genes <- stalk_cell_genes$Feature

# nDCG Calculation
ndcg <- function(relevance, k) {
  # Input: 
  #   relevance - A numeric vector of relevance scores
  #   k - The number of top items to consider
  
  # Truncate relevance to top-k
  relevance <- relevance[1:k]
  
  # Compute DCG
  dcg <- sum(relevance / log2(seq_along(relevance) + 1))
  
  # Compute Ideal DCG (IDCG)
  ideal_relevance <- sort(relevance, decreasing = TRUE)
  idcg <- sum(ideal_relevance / log2(seq_along(ideal_relevance) + 1))
  
  # Calculate nDCG
  ndcg_value <- dcg / idcg
  return(ndcg_value)
}

#this version was 
compute_PR_Chat <- function(graph, genes_of_interest, k, damping ) {
  nodes <- V(graph)$name
  genes_of_interest <- intersect(genes_of_interest, nodes)
  set.seed(139)
  randomized_genes <- sample(genes_of_interest)
  
  folds <- split(randomized_genes, cut(seq_along(randomized_genes), k, labels = FALSE))
  
  aucs <- numeric(k)
  aps <- numeric(k)  # Average Precision
  precisions_at_1000 <- numeric(k)
  recalls_at_1000 <- numeric(k)
  ndcgs <- numeric(k)
  f1_scores<-numeric(k)
  
  for (i in seq_len(k)) {
    test_genes <- folds[[i]]
    training_genes <- setdiff(genes_of_interest, test_genes)
    teleportation_vector <- ifelse(nodes %in% training_genes, 1, 0)
    teleportation_vector <- teleportation_vector / sum(teleportation_vector)
    
    pr <- page_rank(
      graph,
      algo = "prpack",
      directed = FALSE,
      damping = damping,
      personalized = teleportation_vector
    )
    rank <- sort(pr$vector, decreasing = TRUE)
    rank_normalized <- (rank - min(rank)) / (max(rank) - min(rank))
    
    filtered_rank <- rank_normalized[!(names(rank_normalized) %in% training_genes)]
    
    # use the filtered rank from here 
    idx <- which(names(filtered_rank ) %in% test_genes)
    test_labels <- rep(0, length(filtered_rank ))
    test_labels[idx] <- 1
    
    # AUC
    roc_curve_test <- roc(test_labels, filtered_rank)
    aucs[i] <- pROC::auc(roc_curve_test)
    
    # Average Precision (AP)
    pr_curve <- pr.curve(scores.class0 = filtered_rank, weights.class0 = test_labels, curve = TRUE)
    aps[i] <- pr_curve$auc.integral
    
    # Precision@10
    top_k <- names(filtered_rank)[1:1000]
    precisions_at_1000[i] <- length(intersect(top_k, test_genes)) / length(top_k)
    recalls_at_1000[i] <- length(intersect(top_k, test_genes)) / length(test_genes)
    # nDCG
    f1_scores[i] <- 2 * ( precisions_at_1000[i] * recalls_at_1000[i]) / (precisions_at_1000[i] + recalls_at_1000[i])
    
    ########
    ndcgs[i] <- ndcg(test_labels, k = 1000)
  }
  
  #cat(sprintf("Average AUC: %.4f\n", mean(aucs)))
  #cat(sprintf("Average Precision (AP): %.4f\n", mean(aps)))
  #cat(sprintf("Average Precision@100: %.4f\n", mean(precisions_at_1000)))
  #cat(sprintf("Average Precision@100: %.4f\n", mean(recalls_at_1000)))
  #cat(sprintf("Average nDCG@10: %.4f\n", mean(ndcgs)))

  return(data.frame(
    AUC = mean(aucs),
    AP = mean(aps),
    Precision_at_1000 = mean(precisions_at_1000),
    Recalls_at_1000 = mean(recalls_at_1000),
    f1_score = mean(f1_scores),
    nDCG_at_1000 = mean(ndcgs)
  ))
}

compute_PR<-function(graph, genes_of_interest , k){
  
  nodes<-V(graph)$name
  genes_of_interest<-intersect(genes_of_interest , nodes)
  set.seed(139)  # Set a seed for reproducibility
  randomized_genes <- sample(genes_of_interest)
  
  folds <- split(randomized_genes, cut(seq_along(randomized_genes), k, labels = FALSE))
  
  aucs<-c()
  median_ranks<-c()
  average_ranks<-c()
  first_appearances<-c()
  last_appearances<-c()
  for (i in 1:k ){
    test_genes<-folds[[i]]
    training_genes<-setdiff(genes_of_interest, test_genes)
    teleportation_vector <- ifelse(nodes %in% training_genes, 1, 0)
    teleportation_vector <- teleportation_vector / sum(teleportation_vector)
    pr <- page_rank(
      graph,
      algo = "prpack",          # Algorithm used
      directed = FALSE,         # Undirected graph
      damping = 0.1,           # Damping factor (default is 0.85)
      personalized = teleportation_vector
    )
    rank<-sort(pr$vector, decreasing = TRUE)
    
    rank_normalized <- (rank - min(rank)) / (max(rank) - min(rank))
    
    idx<-which (names(rank) %in% test_genes)
    first_appearance<-min(idx)/length(rank)
    last_appearance<-max(idx)/length(rank)
    median_rank<- median(idx)/length(rank)
    average_rank <- mean(idx)/length(rank)
    median_ranks<-c(median_ranks, median_rank)
    average_ranks<-c(average_ranks , average_rank)
    first_appearances<-c(first_appearances , first_appearance)
    last_appearances<-c(last_appearances , last_appearance)
    
    #print(paste("fold " , i , " median rank: " , median_rank , " mean rank: " ,average_rank , " first appearance: " , first_appearance , " laste appearance: " , last_appearance))
    test_labels<-rep(0 , length(rank))
    test_labels[idx]<-1
    roc_curve_test <- roc( test_labels, rank_normalized)
    auc_value_test <- auc(roc_curve_test)
    
    aucs<-c(aucs, auc_value_test)
    #print(paste("test AUC: " , auc_value_test))
  }
  
  df<-data.frame(mean(aucs) , mean(median_ranks) , mean(average_ranks) , mean(first_appearances) , mean(last_appearances))
  colnames(df)<-c("AUC_Value" , "Median_Rank" , "Mean_Rank" , "Mean_First_Appearance" , "Mean_Last_Appearance")
  return(df)
}


run_analysis<-function(interest , threshold , n_folds , write_graph2file = FALSE , damping){
  df_ <- data.frame()
  genes_of_interest <- stalk_cell_genes
  
  for (file in files){
    #print(file)
    study<-unlist(strsplit(file, "_"))[1]
    counts <- read.csv(paste0(datasets_dir , file) )  
    if (dim(counts)[2] < 8) {
      print(paste("skipping the study " , study))
      next
    }
    #print(dim(counts))
    counts <- counts[!duplicated(counts$Feature), ] #removing duplicates
    counts <- na.omit(counts) #removing NAs
    rownames(counts)<-counts$Feature
    counts <- counts[, !(colnames(counts) %in% c("Feature"))]
    if("X" %in% colnames(counts)){
      print(paste("we have X in the file " , file))
    }
    #print(paste(study , dim(counts)) )
    orig_genes<-nrow(counts)

    ######### filtering for expression values and variance #####
    
  ge_threshold<-quantile(rowMeans(counts), probs = seq(0.05, 0.95, by=0.05))[1]
    vars<-apply(counts , 1, var)
    var_threshold <- quantile(vars, probs = seq(0.05, 0.95, by=0.05))[1]
    
    good <- which(
      rowMeans(counts) > ge_threshold &  # Mean expression threshold
      vars > var_threshold  # Variance threshold
    )
    counts <- counts[good, ]
    #print(paste(study , " kept " , length(good)*100.0/orig_genes , "% genes" ))
    ############################################################
    counts <- NormalizeData(counts)
    counts <- na.omit(counts)
    ################################################################
    n_genes<-dim(counts)[1]
    n_samples<-dim(counts)[2]
    # Computing the correlation matrix
    correlation_matrix <- cor(t(counts),
                              method = "pearson")
    
    # Turning NAs into 0
    correlation_matrix[is.na(correlation_matrix)] <- 0
    
    all(rownames(correlation_matrix) == colnames(correlation_matrix))
    
    # Computing the absolute values of the correlation matrix
    correlation_matrix <- abs(correlation_matrix)
    
    diag(correlation_matrix) <- 0  # Remove self-loops
    
    # Create an adjacency matrix based on the threshold
    correlation_matrix <- correlation_matrix > threshold
    
    # Create an undirected graph
    graph <- graph_from_adjacency_matrix(correlation_matrix, mode = "undirected")
    
    #### save graph for further use 
    if(write_graph2file){
      #print(paste("writing fraph for the study " , study ))
      write_graph(graph, file = paste0( graph_outdir , "/", study, "_", threshold,"_graph.graphml" ), format = "graphml")  
    }
    n_nodes<-length(unlist(V(graph)$name))
    n_edges<-length(E(graph))
    #print(paste("n nodes: " , n_nodes ))
    #print(paste("n edges: " , n_edges)) 
    degrees <- degree(graph)
    avg_degree <- mean(degrees)
    #(graph, genes_of_interest, k, damping = 0.85)
    df<-compute_PR_Chat(graph, genes_of_interest,n_folds, damping )

    roww<-data.frame( study , n_genes, n_samples, n_nodes, n_edges, avg_degree , df$AUC , df$AP , df$Precision_at_1000, df$Recalls_at_1000 , df$nDCG_at_1000, df$f1_score , threshold , n_folds , damping)
    colnames(roww)<-c("Dataset","N_genes" , "N_samples", "N_nodes", "N_edges" , "Avg_degree" , "AUC" , "Prec", "Precision_at_1000" , "Recal_at_1000" , "nDCG_at_1000" , "f1_score" , "Cor_thresh" , "n_folds" , "DF")
    
    df_<-rbind(df_, roww)
    #print(study)
    #print(df_)
  }
  return(df_)
}

damping_factors <- seq(0.1, 0.9 , by = 0.1)
#damping_factors<-c(0.9)
for (damping_factor in damping_factors){
  df_<-run_analysis(interest , corr_threshold , n_folds, write_graph , damping_factor )
  file_name<-paste0(outdir , "/" , interest, "_Integrated_" , corr_threshold, "_PR_" , n_folds, "_dump_", damping_factor,  "_Norm_CV_Chat.csv")
  write.csv(file = file_name , df_)
}