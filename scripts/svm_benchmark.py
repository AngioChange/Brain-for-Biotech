import numpy as np
import pandas as pd
import networkx as nx
from sklearn import svm
from sklearn.model_selection import KFold
import useful_functions 
from datetime import datetime


def compute_metrics(predictions, labels):
    predictions_cls = (predictions > 0.5).astype(int)
    accuracy = np.mean(predictions_cls == labels)
    roc_auc = useful_functions.roc_auc_score(labels, predictions)
    precision_recall = useful_functions.precision_recall_curve(labels, predictions)
    pr_auc = useful_functions.auc(precision_recall[1], precision_recall[0])
    return accuracy, roc_auc, pr_auc

def svm_analysis_entire(data_graph, 
                        positive_genes, 
                        output_folder, 
                        best_C, 
                        best_gamma, 
                        best_kernel):
    negative_genes = [gene for gene in data_graph.index if gene not in positive_genes]
    n_permutations = 5
    performance_metrics = []
    accuracies = []
    roc_aucs = []
    pr_aucs = []
    for i in range(n_permutations):
        np.random.shuffle(positive_genes)
        train_size = int(0.8 * len(positive_genes))
        train_genes_positive = positive_genes[:train_size]
        test_genes_positive = positive_genes[train_size:]
        np.random.shuffle(negative_genes)
        train_size = int(0.8 * len(negative_genes))
        train_genes_negative = negative_genes[:train_size]
        test_genes_negative = negative_genes[train_size:]
        train_genes = train_genes_positive + train_genes_negative
        test_genes = test_genes_positive + test_genes_negative      
        print(f"Permutation {i + 1}, Training genes: {len(train_genes)}, Test genes: {len(test_genes)}")
        labels_train = np.array([1 if gene in positive_genes else 0 for gene in train_genes])
        labels_test = np.array([1 if gene in positive_genes else 0 for gene in test_genes])
        clf = svm.SVC(C = best_C, kernel = best_kernel, gamma = best_gamma, probability=True)
        data_graph_train = data_graph.loc[train_genes, :]
        data_graph_test = data_graph.loc[test_genes, :]
        clf.fit(data_graph_train, labels_train)
        prediction = clf.predict_proba(data_graph_test)[:, 1]
        accuracy, roc_auc, pr_auc = compute_metrics(prediction, labels_test)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
    performance_metrics.append({
        "accuracy": np.mean(accuracies),
        "roc_auc": np.mean(roc_aucs),
        "pr_auc": np.mean(pr_aucs),
        "std_accuracy": np.std(accuracies),
        "std_roc_auc": np.std(roc_aucs),
        "std_pr_auc": np.std(pr_aucs)
    })
    performance_df = pd.DataFrame(performance_metrics)
    performance_file = f"{output_folder}/svm_performance_trained_on_subset.csv"
    performance_df.to_csv(performance_file, index=False, sep = ",")
    # train on the whole data and predict on the whole data
    predictions = np.zeros((data_graph.shape[0], n_permutations))
    all_genes = data_graph.index.tolist()
    all_labels = np.array([1 if gene in positive_genes else 0 for gene in all_genes])
    accuracies = []
    roc_aucs = []
    pr_aucs = []
    for i in range(n_permutations):
        clf = svm.SVC(C = best_C, kernel = best_kernel, gamma = best_gamma, probability=True)
        clf.fit(data_graph, all_labels)
        predictions[:, i] = clf.predict_proba(data_graph)[:, 1]
        accuracy, roc_auc, pr_auc = compute_metrics(predictions[:, i], all_labels)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
    performance_metrics.append({
        "accuracy": np.mean(accuracies),
        "roc_auc": np.mean(roc_aucs),
        "pr_auc": np.mean(pr_aucs),
        "std_accuracy": np.std(accuracies),
        "std_roc_auc": np.std(roc_aucs),
        "std_pr_auc": np.std(pr_aucs)
    })
    performance_df = pd.DataFrame(performance_metrics)
    performance_file = f"{output_folder}/svm_performance_entire_final.csv"
    performance_df.to_csv(performance_file, index=False)
    predictions_df = pd.DataFrame(predictions, index=all_genes)
    predictions_df['average_prediction'] = np.mean(predictions, axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
    predictions_file = f"{output_folder}/svm_predictions_entire.csv"
    predictions_df.to_csv(predictions_file, index=False)
    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_file = f"{output_folder}/svm_predictions_only_entire.csv"
    predictions_only_df.to_csv(predictions_only_file, index=False)
    print(f"Final - Average Accuracy: {np.mean(accuracies)}, Average ROC AUC: {np.mean(roc_aucs)}, Average PR AUC: {np.mean(pr_aucs)}")


def svm_analysis(data_graph, training_genes , output_folder):
    # perform SVM classification on the data_graph with labels
    # use cross-validation to evaluate the model
    date_str = datetime.now().strftime('%Y%m%d')
    labels  = None
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    # randomly select 90% of the training genes for training and 10% for testing
    n_permutations = 5
    performance_metrics = []
    # training genes are positive examples
    all_genes = data_graph.index.tolist()
    negative_examples = [gene for gene in all_genes if gene not in training_genes]
    # Iterate through the KFold splits
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
            # and shuffle them to get random genes
            np.random.shuffle(negative_examples)
            # take the same number of negative examples for training to keep the balance
            negative_examples_train = negative_examples[:len(train_genes)]
            train_genes_p.extend(negative_examples_train)
            # the test genes will be all genes besides the training genes
            # thus that is the fold 5 plus the other gene not used for training
            test_genes_p = [gene for gene in all_genes if gene not in train_genes_p]
            print(F"Training genes: {len(train_genes_p)}, Test genes: {len(test_genes_p)}")
            train_labels = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
            test_labels = np.array([1 if gene in training_genes else 0 for gene in test_genes_p])
            # fit the SVM model on the training genes
            data_graph_train = data_graph.loc[train_genes_p, :]
            data_graph_test = data_graph.loc[test_genes_p, :]

            clf = svm.SVC(kernel='linear', probability=True)
            clf.fit(data_graph_train, train_labels)
            # get probabilities for the positive class of the test set
            prediction_proba = clf.predict_proba(data_graph_test)[:, 1]
            accuracy, roc_auc, pr_auc = compute_metrics(prediction_proba, test_labels)
            accuracies.append(accuracy)
            roc_aucs.append(roc_auc)
            pr_aucs.append(pr_auc)
            # print the performance metrics for this permutation and fold
            print(f"Permutation {_ + 1}, Fold {fold_idx + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
        # Store the average performance metrics for this fold
        avg_accuracy = np.mean(accuracies)
        avg_roc_auc = np.mean(roc_aucs)
        avg_pr_auc = np.mean(pr_aucs)

        print(f"Fold {fold_idx + 1} - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
        # Append the performance metrics for this
        performance_metrics.append({
            'fold': fold_idx + 1,  # Use fold_idx + 1 for 1-based fold number
            'accuracy': np.mean(accuracies),
            'roc_auc': np.mean(roc_aucs),
            'pr_auc': np.mean(pr_aucs), 
            'std_accuracy': np.std(accuracies),
            'std_roc_auc': np.std(roc_aucs),
            'std_pr_auc': np.std(pr_aucs)
        })
    
    # convert performance metrics to a DataFrame
    performance_df = pd.DataFrame(performance_metrics)
    # Save the performance metrics to a CSV file with date
    
    performance_file = f"{output_folder}svm_performance_{date_str}.csv"
    performance_df.to_csv(performance_file, index=False)
    print(f"Performance metrics saved to {performance_file}")

    all_labels = np.array([1 if gene in training_genes else 0 for gene in all_genes])

    accuracies = []
    roc_aucs = []
    pr_aucs = []
    predictions = np.zeros((len(all_genes), n_permutations))  # Initialize predictions array
    for i in range(n_permutations):
        np.random.shuffle(negative_examples)
        other_genes_train = negative_examples[:len(training_genes)]
        train_genes_p = training_genes + other_genes_train
        labels = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
        clf = svm.SVC(kernel='linear', probability = True)
        data_graph_p = data_graph.loc[train_genes_p, :]
        clf.fit(data_graph_p, labels)

        prediction = clf.predict_proba(data_graph)[:, 1]  # Get probabilities for the positive class
         # Store the predictions for this permutation
        predictions[:, i] = prediction
        accuracy, roc_auc, pr_auc = compute_metrics(prediction, all_labels)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
       

    predictions_df = pd.DataFrame(predictions, index=all_genes, columns=[f'perm_{i}' for i in range(n_permutations)])

    predictions_df['average_prediction'] = predictions.mean(axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    # order according to the true label in descending order
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
  
    predictions_file = f"{output_folder}svm_predictions_{date_str}.csv"
    predictions_df.to_csv(predictions_file)

    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    #sort the predictions only by the average prediction
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_df.to_csv(f"{output_folder}svm_predictions_only_{date_str}.csv")
    print(f"Predictions saved to {predictions_file}")
    performance_summary = {
        'average_accuracy': np.mean(accuracies),
        'average_roc_auc': np.mean(roc_aucs),
        'average_pr_auc': np.mean(pr_aucs),
        'std_accuracy': np.std(accuracies),
        'std_roc_auc': np.std(roc_aucs),
        'std_pr_auc': np.std(pr_aucs)
    }
    performance_summary_df = pd.DataFrame([performance_summary])
    performance_summary_file = f"{output_folder}svm_performance_summary_{date_str}.csv"
    performance_summary_df.to_csv(performance_summary_file, index=False)
    print(f"Performance summary saved to {performance_summary_file}")
    # Save the SVM model to a file

def run_svm(data_graph, positive_genes, output_folder, entire):
    if entire == "yes":
        svm_analysis_entire(data_graph, positive_genes, output_folder)
    else:
        svm_analysis(data_graph, positive_genes, output_folder)


def svm_5_fold_cv(data_graph, valid_genes, output_folder, best_C, best_gamma, best_kernel):
    # perform SVM classification on the data_graph with labels
    # use cross-validation to evaluate the model
    labels  = None
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    # randomly select 90% of the training genes for training and 10% for testing
    n_permutations = 5
    performance_metrics = []
    # training genes are positive examples
    all_genes = data_graph.index.tolist()
    negative_examples = [gene for gene in all_genes if gene not in training_genes]
    # Iterate through the KFold splits
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
            # and shuffle them to get random genes
            np.random.shuffle(negative_examples)
            # take the same number of negative examples for training to keep the balance
            negative_examples_train = negative_examples[:len(train_genes)]
            train_genes_p.extend(negative_examples_train)
            # the test genes will be all genes besides the training genes
            # thus that is the fold 5 plus the other gene not used for training
            test_genes_p = [gene for gene in all_genes if gene not in train_genes_p]
            print(F"Training genes: {len(train_genes_p)}, Test genes: {len(test_genes_p)}")
            train_labels = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
            test_labels = np.array([1 if gene in training_genes else 0 for gene in test_genes_p])
            # fit the SVM model on the training genes
            data_graph_train = data_graph.loc[train_genes_p, :]
            data_graph_test = data_graph.loc[test_genes_p, :]

            clf = svm.SVC(C = best_C, kernel = best_kernel, gamma = best_gamma, probability=True)
            clf.fit(data_graph_train, train_labels)
            # get probabilities for the positive class of the test set
            prediction_proba = clf.predict_proba(data_graph_test)[:, 1]
            accuracy, roc_auc, pr_auc = compute_metrics(prediction_proba, test_labels)
            accuracies.append(accuracy)
            roc_aucs.append(roc_auc)
            pr_aucs.append(pr_auc)
            # print the performance metrics for this permutation and fold
            print(f"Permutation {_ + 1}, Fold {fold_idx + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
        # Store the average performance metrics for this fold
        avg_accuracy = np.mean(accuracies)
        avg_roc_auc = np.mean(roc_aucs)
        avg_pr_auc = np.mean(pr_aucs)

        print(f"Fold {fold_idx + 1} - Average Accuracy: {avg_accuracy}, Average ROC AUC: {avg_roc_auc}, Average PR AUC: {avg_pr_auc}")
        # Append the performance metrics for this
        performance_metrics.append({
            'fold': fold_idx + 1,  # Use fold_idx + 1 for 1-based fold number
            'accuracy': np.mean(accuracies),
            'roc_auc': np.mean(roc_aucs),
            'pr_auc': np.mean(pr_aucs), 
            'std_accuracy': np.std(accuracies),
            'std_roc_auc': np.std(roc_aucs),
            'std_pr_auc': np.std(pr_aucs)
        })
    
    # convert performance metrics to a DataFrame
    performance_df = pd.DataFrame(performance_metrics)
    # Save the performance metrics to a CSV file with date
    
    performance_file = f"{output_folder}/svm_performance.csv"
    performance_df.to_csv(performance_file, index=False)
    print(f"Performance metrics saved to {performance_file}")

    all_labels = np.array([1 if gene in training_genes else 0 for gene in all_genes])

    accuracies = []
    roc_aucs = []
    pr_aucs = []
    predictions = np.zeros((len(all_genes), n_permutations))  # Initialize predictions array
    for i in range(n_permutations):
        np.random.shuffle(negative_examples)
        other_genes_train = negative_examples[:len(training_genes)]
        train_genes_p = training_genes + other_genes_train
        labels = np.array([1 if gene in training_genes else 0 for gene in train_genes_p])
        clf = svm.SVC(C = best_C, kernel = best_kernel, gamma = best_gamma, probability=True)
        data_graph_p = data_graph.loc[train_genes_p, :]
        clf.fit(data_graph_p, labels)

        prediction = clf.predict_proba(data_graph)[:, 1]  # Get probabilities for the positive class
         # Store the predictions for this permutation
        predictions[:, i] = prediction
        accuracy, roc_auc, pr_auc = compute_metrics(prediction, all_labels)
        accuracies.append(accuracy)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        print(f"Permutation {i + 1}, Accuracy: {accuracy}, ROC AUC: {roc_auc}, PR AUC: {pr_auc}")
       

    predictions_df = pd.DataFrame(predictions, index=all_genes, columns=[f'perm_{i}' for i in range(n_permutations)])

    predictions_df['average_prediction'] = predictions.mean(axis=1)
    predictions_df['true_label'] = all_labels
    predictions_df['gene'] = predictions_df.index
    # order according to the true label in descending order
    predictions_df = predictions_df.sort_values(by='true_label', ascending=False)
    predictions_df = predictions_df.reset_index(drop=True)
  
    predictions_file = f"{output_folder}/svm_predictions.csv"
    predictions_df.to_csv(predictions_file)

    predictions_only_df = predictions_df[(predictions_df['true_label'] == 0) & (predictions_df['average_prediction'] > 0.5)].drop(columns=['true_label'])
    #sort the predictions only by the average prediction
    predictions_only_df = predictions_only_df.sort_values(by='average_prediction', ascending=False)
    predictions_only_df = predictions_only_df.reset_index(drop=True)
    predictions_only_df.to_csv(f"{output_folder}/svm_predictions_only.csv")
    print(f"Predictions saved to {predictions_file}")
    performance_summary = {
        'average_accuracy': np.mean(accuracies),
        'average_roc_auc': np.mean(roc_aucs),
        'average_pr_auc': np.mean(pr_aucs),
        'std_accuracy': np.std(accuracies),
        'std_roc_auc': np.std(roc_aucs),
        'std_pr_auc': np.std(pr_aucs)
    }
    performance_summary_df = pd.DataFrame([performance_summary])
    performance_summary_file = f"{output_folder}/svm_performance_summary.csv"
    performance_summary_df.to_csv(performance_summary_file, index=False)
    print(f"Performance summary saved to {performance_summary_file}")
    # Save the SVM model to a file