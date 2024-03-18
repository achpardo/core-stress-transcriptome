#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import random
import json
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import RandomizedSearchCV
from imblearn.over_sampling import SMOTE
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import normalize
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance


def load_clean_data(path_to_tpm,single_stress="none"):
    """
    Args:
        path_to_tpm = full path to file containing log transformed, BioProject-corrected TPM, with columns for Sample, BioProject, & Treatment
        single_stress = a single stressor to which the data must be subsetted, or "none" (default)
    """
    # load the TPM data
    raw_tpm = pd.read_csv(path_to_tpm,sep="\t",header="infer")
    # if there is a single stress to subset to, subset to that stress
    if single_stress != "none":
        if single_stress == "Drought":
            raw_tpm = raw_tpm[raw_tpm["Treatment"].isin(["Drought","DroughtRepeat","Control"])]
        else:
            raw_tpm = raw_tpm[raw_tpm["Treatment"].isin([single_stress,"Control"])]
        # remove BioProjects that only have control samples left
        testdf = raw_tpm[["Sample","BioProject","Treatment"]]
        for b in testdf["BioProject"].unique():
            df = testdf[testdf["BioProject"]==b]
            if len(df["Treatment"].unique())==1:
                testdf = testdf[testdf["BioProject"]!=b]
        raw_tpm = raw_tpm.merge(testdf)
    # labeling: set Control to 0 and any stress to 1
    proxy = []
    for i in range(len(raw_tpm.index)):
        if raw_tpm.iloc[i,raw_tpm.columns.get_loc("Treatment")] == "Control":
            proxy.append(0)
        else:
            proxy.append(1)
    raw_tpm["Label"] = proxy
    # return the dataframe
    return raw_tpm

def check_if_balanced(labeled_tpm):
    """
    Args:
        labeled_tpm = raw TPM with columns for Sample, BioProject, Label, Treatment
    """
    if labeled_tpm["Label"].value_counts()[0] == labeled_tpm["Label"].value_counts()[1]:
        return True
    else:
        return False

def downsample(dataframe):
    """
    Args:
        dataframe = a log TPM dataframe with a Label column and Sample set as the index
    """
    # generate a variable of value counts
    vc = dataframe["Label"].value_counts()

    # subset data to only samples labeled 1
    ones_only = dataframe[dataframe["Label"]==1]
    
    # downsample from the subsetted dataframe
    ds = ones_only.sample(n=vc[1],random_state=42)

    # subset original data to control samples
    zeroes = dataframe[dataframe["Label"]==0]

    # concatenate controls and downsampled stress samples
    downsampled = pd.concat([ds,zeroes])
    # return dataframe
    return downsampled

def pre_split_process(log_tpm,balanced,downsample=False):
    """
    Args:
        log_tpm = dataframe containing raw TPM values, columns for Sample, BioProject, Treatment, Label
        balanced = Boolean variable, True or False (result of check_if_balanced())
        downsample = Boolean variable, True or False, default False (set manually outside function)
    """
    # temporarily, set index to Sample and drop BioProject, Label, & Treatment columns
    blt = log_tpm[["Sample","BioProject","Treatment","Label"]]
    tpmi = log_tpm.set_index("Sample").drop(["BioProject","Treatment","Label"],axis=1)
    # downsample data if needed
    if balanced==True:
        if downsample==True:
            # add back labels
            tpmi = blt[["Sample","Label"]].merge(tpmi.reset_index().rename(columns={"index":"Sample"}))
            # set Sample as index
            tpmi = tpmi.set_index("Sample")
            # downsample the data
            tpmi = downsample(tpmi)
    # add treatment, labels, and BioProject back in, set Sample as the index again
    labeled = blt.merge(tpmi.reset_index().rename(columns={"index":"Sample"}))
    labeled.set_index("Sample",inplace=True)
    # return dataframe
    return labeled

def split_prep_bioproject(bioproject,dataframe,balance="up"):
    """
    Args:
        bioproject = BioProject to hold out for testing (PRJNAXXXXXX)
        dataframe = starting dataframe of log TPM with labels
        balance = str: "none","up" (downsampling will be done before splitting, outside of this function)
    """
    # split training and testing sets
    test = dataframe[dataframe["BioProject"]==bioproject]
    train = dataframe[dataframe["BioProject"]!=bioproject]
    # drop BioProject column from both sets
    test = test.drop("BioProject",axis=1)
    train = train.drop("BioProject",axis=1)
    # generate X_train, X_test, y_train, and y_test
    ## where X = gene expression values and y = class labels
    train_X = train.drop("Label",axis=1)
    y_train = train["Label"]
    test_X = test.drop("Label",axis=1)
    y_test = test["Label"]
    # if upsampling: do the upsampling using SMOTE
    if balance=="up":
        sm = SMOTE(random_state=42)
        train_X, y_train = sm.fit_resample(train_X,y_train)
    # for X_train and X_test: scale data to a z-score
    scalar = StandardScaler()
    X_train = scalar.fit_transform(train_X)
    X_test = scalar.fit_transform(test_X)
    # return training and test data
    return X_train, y_train, X_test, y_test

def split_prep_stressor(stressor,dataframe,balance="Up"):
    """
    Args:
        stressor = stressor to hold out for testing (all BioProjects)
        dataframe = log TPM dataframe with Sample, Label, BioProject, Treatment columns (or Sample as index)
        balance = str: "none","up" (downsampling will be done before splitting, outside of this function)
    """
    # in case Sample isn't already a column, reset the index and rename the column to Sample
    if "Sample" not in dataframe.columns:
        dataframe = dataframe.reset_index().rename(columns={"index":"Sample"})
    # generate list of unique BioProjects containing the test stressor
    sbp = dataframe[dataframe["Treatment"]==stressor]["BioProject"].unique()
    # split test from train data
    test = dataframe[dataframe["BioProject"].isin(sbp)]
    test = test[test["Treatment"].isin([stressor,"Control"])]
    # pull out training data
    train = dataframe[~dataframe["Sample"].isin(test["Sample"])]
    # for both sets, make Sample the index again
    test = test.set_index("Sample")
    train = train.set_index("Sample")
    # drop BioProject and Treatment columns from both sets
    test = test.drop(["BioProject","Treatment"],axis=1)
    train = train.drop(["BioProject","Treatment"],axis=1)
    # generate X_train, X_test, y_train, and y_test
    ## where X = gene expression values and y = class labels
    train_X = train.drop("Label",axis=1)
    y_train = train["Label"]
    test_X = test.drop("Label",axis=1)
    y_test = test["Label"]
    # if upsampling: do the upsampling using SMOTE
    if balance=="up":
        sm = SMOTE(random_state=42)
        train_X, y_train = sm.fit_resample(train_X,y_train)
    # for X_train and X_test: scale data to a z-score
    scalar = StandardScaler()
    X_train = scalar.fit_transform(train_X)
    X_test = scalar.fit_transform(test_X)
    # return training and test data
    return X_train, y_train, X_test, y_test

def get_tuned_rf(X_train, y_train, random_grid):
    rf = RandomForestClassifier()
    rf_random = RandomizedSearchCV(estimator=rf,
                                  param_distributions=random_grid,
                                  n_iter=150,
                                  cv=5,
                                  verbose=2,
                                  random_state=42,
                                  n_jobs=-1)
    rf_random.fit(X_train, y_train)
    hyper = rf_random.best_params_
    rfclf_tune = RandomForestClassifier(n_estimators=hyper["n_estimators"],
                                min_samples_split=hyper["min_samples_split"],
                                    min_samples_leaf=hyper["min_samples_leaf"],
                                   max_features=hyper["max_features"],
                                   max_depth=hyper["max_depth"],
                                   bootstrap=hyper["bootstrap"])
    return rfclf_tune

def get_scores(y_test,y_pred,filename,AUC):
    """
    Args:
        y_test = true labels of test data
        y_pred = predicted labels of test data
        filename = file name with extension (JSON)
        AUC = float, AUC score output from ROC curve calculations
    """
    f1 = list(f1_score(y_test,y_pred,average=None))
    prec = list(precision_score(y_test,y_pred,average=None))
    rec = list(recall_score(y_test,y_pred,average=None))
    # construct dictionary of F1 and accuracy scores, precision, and recall
    scores = {"Accuracy":accuracy_score(y_test,y_pred),"F1_class_0":f1[0],"F1_class_1":f1[1],
            "Precision_class_0":prec[0],"Precision_class_1":prec[1],"Recall_class_0":rec[0],"Recall_class_1":rec[1],"AUC":AUC}
    # write dictionary to a JSON file
    with open(filename,"w+") as outfile:
        json.dump(scores,outfile,indent=4)

def plot_roc_curve(rf_tuned_fit,X_test,y_test,plot_filename,data_filename):
    """
    Args:
        rf_tuned_fit = tuned and fit random forest model
        X_test = gene expression test data
        y_test = labels for test data
        plot_filename = filename with directory (no extension) for plot
        data_filename = filename with directory and extension (tab delimited) for false & true positive rate data
    """
    # Calculate the probabilities of the classes
    y_prob = rf_tuned_fit.predict_proba(X_test)[:, 1]

    # Compute the ROC curve
    fpr, tpr, _ = roc_curve(y_test,y_prob)
    roc_auc = roc_auc_score(y_test,y_prob)

    # Save ROC curve data: false and true positive rates
    pd.DataFrame([fpr,tpr]).transpose().rename(columns={0:"False Positive Rate",1:"True Positive Rate"}).to_csv(data_filename,sep="\t",header=True,index=False)

    # Plot the ROC curve
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.savefig(plot_filename+".png")
    plt.savefig(plot_filename+".pdf")

    # return AUC score for later use
    return roc_auc

def plot_confusion_matrix(cm,filename,classes=[0,1],normalize=False,cmap=plt.cm.Blues):
    """
    Args:
        cm = confusion matrix
        classes = a list of the classes, default [0,1]
        normalize = Boolean variable stating whether to normalize, default False
        cmap = colormap to use, default Blues
        filename = filename with directory path, no extension
    """
    title = 'Confusion matrix'
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        title = 'Normalized ' + title
    
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in np.ndindex(cm.shape):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.savefig(filename+".png")
    plt.savefig(filename+".pdf")

def extract_importances(rf_tuned_fit,log_tpm,filename):
    """
    Args:
       rf_tuned_fit = tuned and fit random forest model
       log_tpm = log transformed and scaled TPM dataframe
       filename = path to output file, with extension (tab delimited)
    """
    # pull out gene expression values
    X = log_tpm.drop(["BioProject","Treatment","Label"],axis=1)

    # Get feature importances
    feature_importances = rf_tuned_fit.feature_importances_
    # Create a list of (feature, importance) tuples
    feature_importance_list = list(zip(X.columns, feature_importances))
    # Sort the list by importance in descending order
    sorted_feature_importance_list = sorted(feature_importance_list, key=lambda x: x[1], reverse=True)

    # generate dataframe and save file
    pd.DataFrame(sorted_feature_importance_list,columns=["GeneID","Feature_Importance"]).to_csv(filename,sep="\t",header=True,index=False)

def generate_basename(stressor,sampling):
    """
    Args:
        stressor = stressor to be held out for testing
        sampling = Up, Down, or none
    """
    basename = stressor+"Test_"+sampling
    return basename

def main():
    parser = argparse.ArgumentParser(description="Parse args")
    parser.add_argument("--tpm_file","-t",type=str,help="full path to TPM file")
    parser.add_argument("--single_stress","-i",type=str,help="options: none or a single stress, e.g. Drought")
    parser.add_argument("--sampling","-s",type=str,help="whether to upsample, downsample, or do nothing if data are unbalanced, options: Up, Down, or none")
    parser.add_argument("--path","-p",type=str,help="path to parent directory where output directory should be placed")
    args = parser.parse_args()
    tpm_file = str(args.tpm_file)
    single_stress = str(args.single_stress)
    sampling = str(args.sampling)
    dirpath = str(args.path)

    # load and clean the TPM data; subset to a single stress if specified
    cleaned_tpm = load_clean_data(tpm_file,single_stress="none")

    # set variables for downsampling and whether data are balanced
    bal = check_if_balanced(cleaned_tpm)
    if sampling == "Down":
        ds = True
        us = "none"
    elif sampling == "Up":
        ds = False
        us = "Up"
    else:
        ds = False
        us = "none"

    # log transform the data, and downsample if specified
    log_tpm = pre_split_process(cleaned_tpm,bal,ds)

    # split the data into training and testing sets
    X_train, y_train, X_test, y_test = split_prep_stressor(single_stress,log_tpm,us)

    # create dictionary of hyperparameters to search
    random_search_grid = {'bootstrap': [True, False],
        'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
        'max_features': ['sqrt','log2','none'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]}

    # tune the model
    rfclf_tuned = get_tuned_rf(X_train,y_train,random_search_grid)

    # fit the model again on the training data
    rfclf_tuned.fit(X_train,y_train)

    # make predictions on testing set
    y_pred = rfclf_tuned.predict(X_test)

    # Generate the output directory
    base = generate_basename(single_stress,sampling)
    directory = dirpath+"/"+base+"_output"
    if not os.path.exists(directory):
        os.mkdir(directory)

    # output ROC curve
    auc = plot_roc_curve(rfclf_tuned,X_test,y_test,directory+"/"+base+"_ROC_Curve",directory+"/"+base+"_FPR_TPR_data.tsv")
    plt.clf()

    # output accuracy, F1, precision, and recall as a file
    get_scores(y_test,y_pred,directory+"/"+base+"_scores.json",auc)

    # output confusion matrix
    cm = confusion_matrix(y_test,y_pred)
    plot_confusion_matrix(cm,directory+"/"+base+"_ConfusionMatrix")

    # output feature importances
    extract_importances(rfclf_tuned,log_tpm,directory+"/"+base+"_Feature_Importances_sorted.tsv")

if __name__ == "__main__":
    main()
