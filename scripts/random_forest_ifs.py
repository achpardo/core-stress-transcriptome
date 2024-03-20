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

def subset_tpm(raw_tpm,gene_list):
    # pull out all columns of raw_tpm that do not start with Zm
    cols = []
    for c in raw_tpm.columns:
        if not c.startswith("Zm"):
            cols.append(c)
    glist = gene_list + cols
    subtpm = raw_tpm[glist]
    return subtpm

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

def main():
    parser = argparse.ArgumentParser(description="Parse args")
    parser.add_argument("--tpm_file","-t",type=str,help="full path to TPM file")
    parser.add_argument("--single_stress","-i",type=str,help="options: none or a single stress, e.g. Drought")
    parser.add_argument("--sampling","-s",type=str,help="whether to upsample, downsample, or do nothing if data are unbalanced, options: Up, Down, or none")
    parser.add_argument("--path","-p",type=str,help="path to parent directory where output directory should be placed")
    parser.add_argument("--features_json","-f",type=str,help="path to JSON with subsets of features")
    parser.add_argument("--dataset","-d",type=str,help="options: All or Leaf")
    args = parser.parse_args()
    tpm_file = str(args.tpm_file)
    single_stress = str(args.single_stress)
    sampling = str(args.sampling)
    dirpath = str(args.path)
    features_json = str(args.features_json)

    # load and clean the TPM data; subset to a single stress if specified
    cleaned_tpm = load_clean_data(tpm_file,single_stress="none")

    # load the feature list json file
    fjson = json.load(open(features_json))
    # subset TPM to certain features only; create a dictionary of output dataframes
    subtpmdict = {}
    for k,v in fjson.items():
        print(k)
        subtpmdict[k] = subset_tpm(cleaned_tpm,v)

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

    # process the data, and downsample if specified, for each entry in subtpmdict
    log_subtpmdict = {}
    for k,v in subtpmdict.items():
        log_tpm = pre_split_process(v,bal,ds)
        log_subtpmdict[k] = log_tpm

    # create dictionary of hyperparameters to search
    random_search_grid = {'bootstrap': [True, False],
        'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
        'max_features': ['sqrt','log2','none'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]}

    # loop through log_subtpmdict to run random forest
    ## outputs to save: accuracy, AUC, F1
    outputsdict = {"n_features":[],"accuracy":[],"AUC":[],"F1_class0":[],"F1_class1":[]}
    for k,v in log_subtpmdict.items():
        # split the data into training and testing sets
        X_train, y_train, X_test, y_test = split_prep_stressor(single_stress,v,us)

        # tune the model
        rfclf_tuned = get_tuned_rf(X_train,y_train,random_search_grid)

        # fit the model again on the training data
        rfclf_tuned.fit(X_train,y_train)

        # make predictions on testing set
        y_pred = rfclf_tuned.predict(X_test)
        
        # save number of features
        outputsdict["n_features"].append(k)

        # get & save accuracy
        acc = accuracy_score(y_test,y_pred)
        outputsdict["accuracy"] = acc
        
        # get & save F1 scores
        f1 = list(f1_score(y_test,y_pred,average=None))
        outputsdict["F1_class0"] = f1[0]
        outputsdict["F1_class1"] = f1[1]

        # get & save AUC
        y_prob = rfclf_tuned.predict_proba(X_test)[:,1]
        roc_auc = roc_auc_score(y_test,y_prob)
        outputsdict["AUC"] = roc_auc

    # convert outputsdict to a dataframe
    outdf = pd.DataFrame(outputsdict)

    # save output
    outdf.to_csv(os.path.join(dirpath,single_stress+"_"+dataset+"_IFS_Model_Scores.tsv"),sep="\t",header=True,index=False)

if __name__ == "__main__":
    main()
