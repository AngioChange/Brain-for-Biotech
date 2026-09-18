import numpy as np
import os
import pandas as pd
import networkx as nx
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
import useful_functions    
from datetime import date, datetime
from tqdm import tqdm

# take care, here the predictions are probabilities for the positive class
def compute_metrics(predictions, labels):
    # predictions are probabilites for the positive class
    # labels are the true labels
    predictions_cls = (predictions > 0.5).astype(int)  # Convert probabilities to class labels
    accuracy = np.mean(predictions_cls == labels)
    roc_auc = useful_functions.roc_auc_score(labels, predictions)
    precision_recall = useful_functions.precision_recall_curve(labels, predictions) 
    pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
    return accuracy, roc_auc, pr_auc
    
def random_forest_analysis_entire(data_graph, positive_genes, output_folder, best_n_estimators, best_features, best_depth):
    # perfrom random forest on 4/5 training genes and 1.5 testing genes
    # select 80% genes for training and 20% for testing     
    date_str = datetime.now().strftime('%Y%m%d')
    negative_genes = [gene for gene in data_graph.index if gene not in positive_genes]
    # permute training genes
    
    n_permutations = 5  
    performance_metrics = []
    accuracies = []
    roc_aucs = []
    pr_aucs = []
   
    for i in range(n_permutations):
        # randomly select 80% of the positive genes for training
        np.random.shuffle(positive_genes)
        train_size = int(0.8 * len(positive_genes))
        train_genes_positive  = positive_genes[:train_size]
        test_genes_positive = positive_genes[train_size:]
        # randomly select 80% of the negative genes for training
        np.random.shuffle(negative_genes)
        train_size = int(0.8 * len(negative_genes))
        train_genes_negative = negative_genes[:train_size]
        test_genes_negative = negative_genes[train_size:]
        # combine the training genes
        train_genes = train_genes_positive + train_genes_negative
        test_genes = test_genes_positive + test_genes_negative
        print(f"Permutation {i + 1}, Training genes: {len(train_genes)}, Test genes: {len(test_genes)}")
        # create labels for training and testing
        labels_train = np.array([1 if gene in positive_genes else 0 for gene in train_genes])
        labels_test = np.array([1 if gene in positive_genes else 0 for gene in test_genes])

        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        data_graph_train = data_graph.loc[train_genes, :]
        data_graph_test = data_graph.loc[test_genes, :]
        clf.fit(data_graph_train, labels_train) 
        # no predictions here
        prediction = clf.predict_proba(data_graph_test)[:, 1]  # Get probabilities for the positive class
        accuracy, roc_auc, pr_auc = compute_metrics(prediction, labels_test)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
    # store the performance metrics for all permutations
    performance_metrics.append({
        "accuracy": np.mean(accuracies),
        "roc_auc": np.mean(roc_aucs),
        "pr_auc": np.mean(pr_aucs),
        "std_accuracy": np.std(accuracies),
        "std_roc_auc": np.std(roc_aucs),
        "std_pr_auc": np.std(pr_aucs)
    })
    # convert the performance metrics to a DataFrame
    performance_df = pd.DataFrame(performance_metrics)
    performance_file = f"{output_folder}random_forest_performance_entire_{date_str}.csv"
    performance_df.to_csv(performance_file, index=False)    
    # train on the whole data and predict on the whole data
    predictions = np.zeros((data_graph.shape[0], n_permutations))
    all_genes = data_graph.index.tolist()
    all_labels = np.array([1 if gene in positive_genes else 0 for gene in all_genes])
    accuracies = []
    roc_aucs = []
    pr_aucs = []
    for i in range(n_permutations): 
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(data_graph, all_labels)
        predictions[:, i] = clf.predict_proba(data_graph)[:, 1]
        accuracy, roc_auc, pr_auc = compute_metrics(predictions[:, i], all_labels)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
    # store the performance metrics for the whole data
    performance_metrics.append({
        "accuracy": np.mean(accuracies),
        "roc_auc": np.mean(roc_aucs),
        "pr_auc": np.mean(pr_aucs),
        "std_accuracy": np.std(accuracies),
        "std_roc_auc": np.std(roc_aucs),
        "std_pr_auc": np.std(pr_aucs)
    })
    performance_df = pd.DataFrame(performance_metrics)
    performance_file = f"{output_folder}random_forest_performance_entire_final_{date_str}.csv"
    performance_df.to_csv(performance_file, index=False)
    # save the predictions to a DataFrame
    predictions_df = pd.DataFrame(predictions, index=all_genes)
    predictions_df['average_prediction'] = np.mean(predictions, axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
    predictions_file = f"{output_folder}random_forest_predictions_entire_{date_str}.csv"
    predictions_df.to_csv(predictions_file, index=False)    
    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    # sort the predictions only by the average prediction
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_file = f"{output_folder}random_forest_predictions_only_entire_{date_str}.csv"
    predictions_only_df.to_csv(predictions_only_file, index=False)  
    print(f"Final - Average Accuracy: {np.mean(accuracies)}, Average ROC AUC: {np.mean(roc_aucs)}, Average PR AUC: {np.mean(pr_aucs)}")

    performance_metrics.append({
        'fold': 'final',  # Use 'final' to indicate the final performance
        'accuracy': np.mean(accuracies),
        'roc_auc': np.mean(roc_aucs),
        'pr_auc': np.mean(pr_aucs),
        'std_accuracy': np.std(accuracies),
        'std_roc_auc': np.std(roc_aucs),
        'std_pr_auc': np.std(pr_aucs)
    })
    performance_df = pd.DataFrame(performance_metrics)
    performance_file = f"{output_folder}random_forest_performance_summary_entire_{date_str}.csv"
    performance_df.to_csv(performance_file, index=False)
    print(f"Performance metrics saved to {performance_file}")

def random_forest_analysis(data_graph, training_genes, output_folder):
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    # randomly select 90% of the training genes for training and 10% for testing
    n_permutations = 5
    performance_metrics = []
    all_genes = data_graph.index.tolist()
    negative_examples = [gene for gene in all_genes if gene not in training_genes]
    all_labels = np.array([1 if gene in training_genes else 0 for gene in all_genes])
    date_str = datetime.now().strftime('%Y%m%d')
    print(f"Number of training genes: {len(training_genes)}")
    for fold_idx, (train_index, test_index) in enumerate(kf.split(training_genes)):    
        train_genes = [training_genes[i] for i in train_index]
 
        # add the same number of genes from all genes besides the traning genes and test genes
        # to balance the dataset
        accuracies = []
        roc_aucs = []
        pr_aucs = []
        for _ in range(n_permutations):
            # first copy the training and test genes otherwise they will be modified
            # during the random selection of other genes
            train_genes_p = train_genes.copy()
            np.random.shuffle(negative_examples)
            other_genes_train = negative_examples[:len(train_genes)]
            train_genes_p.extend(other_genes_train)
            test_genes_p = [gene for gene in all_genes if gene not in train_genes_p]
            labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
            labels_test = np.array([1 if gene in training_genes else 0 for gene in test_genes_p])
            # fit the Random Forest model on the training genes
            data_graph_train = data_graph.loc[train_genes_p, :]
            data_graph_test = data_graph.loc[test_genes_p, :]

            clf = RandomForestClassifier(n_estimators=100, random_state=42)
            clf.fit(data_graph_train, labels_train)
            # Evaluate the model on the test set
         
            accuracy = clf.score(data_graph_test, labels_test)
            accuracies.append(accuracy)
            roc_auc = useful_functions.roc_auc_score(labels_test, clf.predict_proba(data_graph_test)[:, 1])
            roc_aucs.append(roc_auc)
            precision_recall = useful_functions.precision_recall_curve(labels_test, clf.predict_proba(data_graph_test)[:, 1])
            pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
            pr_aucs.append(pr_auc)
            print(f"Permutation {_ + 1}, Fold {kf.get_n_splits()}, Accuracy: {accuracy}")
        # Store the average performance metrics for this fold
        avg_accuracy = np.mean(accuracies)
        avg_roc_auc = np.mean(roc_aucs)
        avg_pr_auc = np.mean(pr_aucs)
        print(f"Fold {fold_idx + 1} - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
        # Append the performance metrics for this
        performance_metrics.append({
            'fold': fold_idx + 1,  # Use fold_idx + 1 for 1-based fold number
            'accuracy': avg_accuracy,
            'roc_auc': avg_roc_auc,
            'pr_auc': avg_pr_auc,
            'std_accuracy': np.std(accuracies),
            'std_roc_auc': np.std(roc_aucs),
            'std_pr_auc': np.std(pr_aucs)
        })
    # conver the performance metrics to a DataFrame
    performance_df = pd.DataFrame(performance_metrics)
    filename = f"{output_folder}random_forest_performance_{datetime.now().strftime('%Y%m%d')}.csv"
    performance_df.to_csv(filename, index=False)
    print(f"Performance metrics saved to {filename}")
    # train 10 times on the whole training set and predict on the whole data
    # this is used to get the final performance metrics
    accuracies = []
    roc_aucs = []
    pr_aucs = []
    predictions = np.zeros((data_graph.shape[0], n_permutations))
    for  i in range(n_permutations):     
        np.random.shuffle(negative_examples)
        other_genes_train = negative_examples[:len(training_genes)]
        train_genes_p = training_genes.copy()
        train_genes_p.extend(other_genes_train)
        labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
        data_graph_train = data_graph.loc[train_genes_p, :]
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(data_graph_train, labels_train)
        # use the model to predict on the whole data the labels are alredy defined
        prediction = clf.predict_proba(data_graph)[:, 1]  # Get probabilities for the positive class
        predictions[:, i] = prediction
        accuracy = clf.score(data_graph, all_labels)
        accuracies.append(accuracy)
        roc_auc = useful_functions.roc_auc_score(all_labels, prediction)
        roc_aucs.append(roc_auc)
        precision_recall = useful_functions.precision_recall_curve(all_labels, prediction)
        pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy} , ROC AUC: {roc_auc}, PR AUC: {pr_auc}")

    # Store the average performance metrics for the whole data
    avg_accuracy = np.mean(accuracies)
    avg_roc_auc = np.mean(roc_aucs)
    avg_pr_auc = np.mean(pr_aucs)       
    predictions_df = pd.DataFrame(predictions, index=all_genes)
    predictions_df['average_prediction'] = np.mean(predictions, axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
    filename = f"{output_folder}random_forest_predictions_{date_str}.csv"
    predictions_df.to_csv(filename, index=False)

    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    #sort the predictions only by the average prediction
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_df.to_csv(f"{output_folder}random_forest_predictions_only_{date_str}.csv", index=False)

    print(f"Final - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
    performance_metrics.append({
        'fold': 'final',  # Use 'final' to indicate the final performance
        'accuracy': avg_accuracy,
        'roc_auc': avg_roc_auc,
        'pr_auc': avg_pr_auc,
        'std_accuracy': np.std(accuracies),
        'std_roc_auc': np.std(roc_aucs),
        'std_pr_auc': np.std(pr_aucs)
    })

    filename = f"{output_folder}random_forest_performance_summary_{date_str}.csv"
    performance_df = pd.DataFrame(performance_metrics)
    performance_df.to_csv(filename, index=False)
    print(f"Performance metrics saved to {filename}")

def run_rf(data_graph, training_genes, output_folder, entire):
    if entire == "yes":
        print("Running Random Forest analysis on the entire dataset")
        random_forest_anaysis_entire(data_graph, training_genes, output_folder)
    else:
        print("Running Random Forest analysis with K-Fold cross-validation")
        random_forest_analysis(data_graph, training_genes, output_folder)


def benchmark_hyperparameters_rf(data_graph, training_genes):
    """
    Computes a benchmark on 3 hyperparameters:
    1) n_estimators (default: 100) [0 - 500 , 100]
    2) max_features (default: "sqrt") ["log2", None]
    3) max_depth (default: None) [0 - 500, 100]
    """
    for n_estimators in tqdm(list(range(1, 501, 100)), desc = "Processing n_estimators ...", position = 0):
        for max_features in tqdm(["log2", None], desc = "Processing max_features ...", position = 1):
            for max_depth in tqdm(list(range(1, 501, 100)), desc = "Processing max_depth", position = 2):
                print(f"N_estimators considered: {n_estimators}")
                print(f"max_features considered: {max_features}")
                print(f"max_depth considered: {max_depth}")

                output_folder = f"/Users/u0148349/Desktop/Miscellaneous/BfBio/RF/Benchmark/N_estimators_{n_estimators}_max_features_{max_features}_max_depth_{max_depth}/"
                os.mkdir(output_folder)
        
                kf = KFold(n_splits=5, shuffle=True, random_state=42)
                # randomly select 90% of the training genes for training and 10% for testing
                n_permutations = 5
                performance_metrics = []
                all_genes = data_graph.index.tolist()
                negative_examples = [gene for gene in all_genes if gene not in training_genes]
                all_labels = np.array([1 if gene in training_genes else 0 for gene in all_genes])
                print(f"Number of training genes: {len(training_genes)}")
                for fold_idx, (train_index, test_index) in enumerate(kf.split(training_genes)):    
                    train_genes = [training_genes[i] for i in train_index]
 
                    # add the same number of genes from all genes besides the traning genes and test genes
                    # to balance the dataset
                    accuracies = []
                    roc_aucs = []
                    pr_aucs = []
                    for _ in range(n_permutations):
                
                        # first copy the training and test genes otherwise they will be modified
                        # during the random selection of other genes
                        train_genes_p = train_genes.copy()
                        np.random.shuffle(negative_examples)
                        other_genes_train = negative_examples[:len(train_genes)]
                        train_genes_p.extend(other_genes_train)
                        test_genes_p = [gene for gene in all_genes if gene not in train_genes_p]
                        labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
                        labels_test = np.array([1 if gene in training_genes else 0 for gene in test_genes_p])
                    
                        # fit the Random Forest model on the training genes
                        data_graph_train = data_graph.loc[train_genes_p, :]
                        data_graph_test = data_graph.loc[test_genes_p, :]

                        clf = RandomForestClassifier(n_estimators=n_estimators,
                                                     max_depth = max_depth,
                                                     max_features = max_features,
                                                     random_state=42)
                
                        clf.fit(data_graph_train, labels_train)
                
                        # Evaluate the model on the test set
                        accuracy = clf.score(data_graph_test, labels_test)
                        accuracies.append(accuracy)
                        roc_auc = useful_functions.roc_auc_score(labels_test, clf.predict_proba(data_graph_test)[:, 1])
                        roc_aucs.append(roc_auc)
                        precision_recall = useful_functions.precision_recall_curve(labels_test,
                                                                               clf.predict_proba(data_graph_test)[:, 1])
                        pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
                        pr_aucs.append(pr_auc)
                        print(f"Permutation {_ + 1}, Fold {kf.get_n_splits()}, Accuracy: {accuracy}")
                
                    # Store the average performance metrics for this fold
                    avg_accuracy = np.mean(accuracies)
                    avg_roc_auc = np.mean(roc_aucs)
                    avg_pr_auc = np.mean(pr_aucs)
                    print(f"Fold {fold_idx + 1} - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
            
                    # Append the performance metrics for this
                    performance_metrics.append({
                        'fold': fold_idx + 1,  # Use fold_idx + 1 for 1-based fold number
                        'accuracy': avg_accuracy,
                        'roc_auc': avg_roc_auc,
                        'pr_auc': avg_pr_auc,
                        'std_accuracy': np.std(accuracies),
                        'std_roc_auc': np.std(roc_aucs),
                        'std_pr_auc': np.std(pr_aucs)
                    })
            
                    # convert the performance metrics to a DataFrame
                    performance_df = pd.DataFrame(performance_metrics)
                    filename = f"{output_folder}random_forest_performance.csv"
                    performance_df.to_csv(filename, index=False)
                    print(f"Performance metrics saved to {filename}")

                # train 10 times on the whole training set and predict on the whole data
                # this is used to get the final performance metrics

                accuracies = []
                roc_aucs = []
                pr_aucs = []
                predictions = np.zeros((data_graph.shape[0], n_permutations))
                for i in range(n_permutations):     
                    np.random.shuffle(negative_examples)
                    other_genes_train = negative_examples[:len(training_genes)]
                    train_genes_p = training_genes.copy()
                    train_genes_p.extend(other_genes_train)
                    labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
                    data_graph_train = data_graph.loc[train_genes_p, :]
                    clf = RandomForestClassifier(n_estimators=n_estimators, 
                                                 max_depth = max_depth,
                                                 max_features = max_features,
                                                 random_state=42)
                
                    clf.fit(data_graph_train, labels_train)
            
                    # use the model to predict on the whole data the labels are alredy defined
                    prediction = clf.predict_proba(data_graph)[:, 1]  # Get probabilities for the positive class
                    predictions[:, i] = prediction
                    accuracy = clf.score(data_graph, all_labels)
                    accuracies.append(accuracy)
                    roc_auc = useful_functions.roc_auc_score(all_labels, prediction)
                    roc_aucs.append(roc_auc)
                    precision_recall = useful_functions.precision_recall_curve(all_labels, prediction)
                    pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
                    pr_aucs.append(pr_auc)
                    print(f"Permutation {i + 1}, Accuracy: {accuracy} , ROC AUC: {roc_auc}, PR AUC: {pr_auc}")

                # Store the average performance metrics for the whole data
                avg_accuracy = np.mean(accuracies)
                avg_roc_auc = np.mean(roc_aucs)
                avg_pr_auc = np.mean(pr_aucs)       
                predictions_df = pd.DataFrame(predictions, index=all_genes)
                predictions_df['average_prediction'] = np.mean(predictions, axis=1)
                predictions_df['true_label'] = all_labels
                predictions_df['gene'] = predictions_df.index
                predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
                predictions_df = predictions_df.reset_index(drop=True)
                filename = f"{output_folder}random_forest_predictions.csv"
                predictions_df.to_csv(filename, index=False)

                predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
                #sort the predictions only by the average prediction
                predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
                predictions_only_df = predictions_only_df.reset_index(drop=True)
                predictions_only_df.to_csv(f"{output_folder}random_forest_predictions_only.csv", index=False)

                print(f"Final - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
                performance_metrics.append({
                    'fold': 'final',  # Use 'final' to indicate the final performance
                    'accuracy': avg_accuracy,
                    'roc_auc': avg_roc_auc,
                    'pr_auc': avg_pr_auc,
                    'std_accuracy': np.std(accuracies),
                    'std_roc_auc': np.std(roc_aucs),
                    'std_pr_auc': np.std(pr_aucs)
                })

                filename = f"{output_folder}random_forest_performance_summary.csv"
                performance_df = pd.DataFrame(performance_metrics)
                performance_df.to_csv(filename, index=False)
                print(f"Performance metrics saved to {filename}")


def random_forest_5Fold_CV(data_graph, training_genes, output_folder, best_n_estimators, best_features, best_depth):
    
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    # randomly select 90% of the training genes for training and 10% for testing
    n_permutations = 5
    performance_metrics = []
    all_genes = data_graph.index.tolist()
    negative_examples = [gene for gene in all_genes if gene not in training_genes]
    all_labels = np.array([1 if gene in training_genes else 0 for gene in all_genes])
    print(f"Number of training genes: {len(training_genes)}")
    for fold_idx, (train_index, test_index) in enumerate(kf.split(training_genes)):    
        train_genes = [training_genes[i] for i in train_index]
 
        # add the same number of genes from all genes besides the traning genes and test genes
        # to balance the dataset
        accuracies = []
        roc_aucs = []
        pr_aucs = []
        for _ in range(n_permutations):
            # first copy the training and test genes otherwise they will be modified
            # during the random selection of other genes
            train_genes_p = train_genes.copy()
            np.random.shuffle(negative_examples)
            other_genes_train = negative_examples[:len(train_genes)]
            train_genes_p.extend(other_genes_train)
            test_genes_p = [gene for gene in all_genes if gene not in train_genes_p]
            labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
            labels_test = np.array([1 if gene in training_genes else 0 for gene in test_genes_p])
            # fit the Random Forest model on the training genes
            data_graph_train = data_graph.loc[train_genes_p, :]
            data_graph_test = data_graph.loc[test_genes_p, :]

            clf = RandomForestClassifier(n_estimators = best_n_estimators,
                                         random_state=42,
                                         max_depth = best_depth,
                                         max_features = best_features)
            
            clf.fit(data_graph_train, labels_train)
            # Evaluate the model on the test set
         
            accuracy = clf.score(data_graph_test, labels_test)
            accuracies.append(accuracy)
            roc_auc = useful_functions.roc_auc_score(labels_test, clf.predict_proba(data_graph_test)[:, 1])
            roc_aucs.append(roc_auc)
            precision_recall = useful_functions.precision_recall_curve(labels_test, clf.predict_proba(data_graph_test)[:, 1])
            pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
            pr_aucs.append(pr_auc)
            print(f"Permutation {_ + 1}, Fold {kf.get_n_splits()}, Accuracy: {accuracy}")
        # Store the average performance metrics for this fold
        avg_accuracy = np.mean(accuracies)
        avg_roc_auc = np.mean(roc_aucs)
        avg_pr_auc = np.mean(pr_aucs)
        print(f"Fold {fold_idx + 1} - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
        # Append the performance metrics for this
        performance_metrics.append({
            'fold': fold_idx + 1,  # Use fold_idx + 1 for 1-based fold number
            'accuracy': avg_accuracy,
            'roc_auc': avg_roc_auc,
            'pr_auc': avg_pr_auc,
            'std_accuracy': np.std(accuracies),
            'std_roc_auc': np.std(roc_aucs),
            'std_pr_auc': np.std(pr_aucs)
        })
    # conver the performance metrics to a DataFrame
    performance_df = pd.DataFrame(performance_metrics)
    filename = f"{output_folder}random_forest_performance.csv"
    performance_df.to_csv(filename, index=False)
    print(f"Performance metrics saved to {filename}")
    # train 10 times on the whole training set and predict on the whole data
    # this is used to get the final performance metrics
    accuracies = []
    roc_aucs = []
    pr_aucs = []
    predictions = np.zeros((data_graph.shape[0], n_permutations))
    for  i in range(n_permutations):     
        np.random.shuffle(negative_examples)
        other_genes_train = negative_examples[:len(training_genes)]
        train_genes_p = training_genes.copy()
        train_genes_p.extend(other_genes_train)
        labels_train = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
        data_graph_train = data_graph.loc[train_genes_p, :]
        clf = RandomForestClassifier(n_estimators = best_n_estimators, 
                                     random_state=42,
                                     max_depth = best_depth,
                                     max_features = best_features)
        
        clf.fit(data_graph_train, labels_train)
        # use the model to predict on the whole data the labels are alredy defined
        prediction = clf.predict_proba(data_graph)[:, 1]  # Get probabilities for the positive class
        predictions[:, i] = prediction
        accuracy = clf.score(data_graph, all_labels)
        accuracies.append(accuracy)
        roc_auc = useful_functions.roc_auc_score(all_labels, prediction)
        roc_aucs.append(roc_auc)
        precision_recall = useful_functions.precision_recall_curve(all_labels, prediction)
        pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy} , ROC AUC: {roc_auc}, PR AUC: {pr_auc}")

    # Store the average performance metrics for the whole data
    avg_accuracy = np.mean(accuracies)
    avg_roc_auc = np.mean(roc_aucs)
    avg_pr_auc = np.mean(pr_aucs)       
    predictions_df = pd.DataFrame(predictions, index=all_genes)
    predictions_df['average_prediction'] = np.mean(predictions, axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
    filename = f"{output_folder}random_forest_predictions.csv"
    predictions_df.to_csv(filename, index=False)

    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    #sort the predictions only by the average prediction
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_df.to_csv(f"{output_folder}random_forest_predictions_only.csv", index=False)

    print(f"Final - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
    performance_metrics.append({
        'fold': 'final',  # Use 'final' to indicate the final performance
        'accuracy': avg_accuracy,
        'roc_auc': avg_roc_auc,
        'pr_auc': avg_pr_auc,
        'std_accuracy': np.std(accuracies),
        'std_roc_auc': np.std(roc_aucs),
        'std_pr_auc': np.std(pr_aucs)
    })

    filename = f"{output_folder}random_forest_performance_summary.csv"
    performance_df = pd.DataFrame(performance_metrics)
    performance_df.to_csv(filename, index=False)
    print(f"Performance metrics saved to {filename}")