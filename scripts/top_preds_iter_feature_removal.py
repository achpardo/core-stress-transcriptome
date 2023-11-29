import pandas as pd
import statistics
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.metrics import classification_report, accuracy_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import confusion_matrix
from imblearn.over_sampling import SMOTE
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import normalize
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.preprocessing import StandardScaler
import json
from sklearn.metrics import f1_score

def load_clean_data(path_to_tpm):
    """
    Args:
        path_to_tpm = full path to file containing raw TPM, columns for Sample, BioProject, & Treatment
        single_stress = a single stressor to which the data must be subsetted, or "none" (default)
    """
    # load the TPM data
    raw_tpm = pd.read_csv(path_to_tpm,sep="\t",header="infer")
    # replace DroughtRepeat with Drought
    raw_tpm["Treatment"].mask(raw_tpm["Treatment"]=="DroughtRepeat","Drought",inplace=True)
    # labeling: set Control to 0 and stressors to 1-6
    slabels = {"Control":0,"Drought":1,"Cold":2,"Heat":3,"Flooding":4,"Low_Nitrogen":5,"Salt":6}
    proxy = []
    for i in range(len(raw_tpm.index)):
        t = raw_tpm.iloc[i,raw_tpm.columns.get_loc("Treatment")]
        proxy.append(slabels[t])
    raw_tpm["Label"] = proxy
    # return the dataframe
    return raw_tpm

def variance_threshold_selector(data):
    selector = VarianceThreshold()
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]

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

def check_if_balanced(labeled_tpm):
    """
    Args:
        labeled_tpm = raw TPM with columns for Sample, BioProject, Label, Treatment
    """
    if labeled_tpm["Label"].value_counts()[0] == labeled_tpm["Label"].value_counts()[1]:
        return True
    else:
        return False

def pre_split_transform(raw_tpm,balanced,downsample=False):
    """
    Args:
        raw_tpm = dataframe containing raw TPM values, columns for Sample, BioProject, Treatment, Label
        balanced = Boolean variable, True or False (result of check_if_balanced())
        downsample = Boolean variable, True or False, default False (set manually outside function)
    """
    # if data have treatment column, drop it
    #if "Treatment" in raw_tpm.columns:
    #    raw_tpm = raw_tpm.drop("Treatment",axis=1)
    # temporarily, set index to Sample and drop BioProject, Label, & Treatment columns
    blt = raw_tpm[["Sample","BioProject","Treatment","Label"]]
    tpmi = raw_tpm.set_index("Sample").drop(["BioProject","Treatment","Label"],axis=1)
    # remove zero-variance genes
    vttpm = variance_threshold_selector(tpmi)
    # log-transform TPM
    vttpm_log = vttpm.apply(lambda x: np.log2(x+1))
    # downsample data if needed
    if balanced!=True:
        if downsample==True:
            # add back labels
            vttpm_log = blt[["Sample","Label"]].merge(vttpm_log.reset_index().rename(columns={"index":"Sample"}))
            # set Sample as index
            vttpm_log = vttpm_log.set_index("Sample")
            # downsample the data
            vttpm_log = downsample(vttpm_log)
    # add treatment, labels, and BioProject back in, set Sample as the index again
    labeled = blt.merge(vttpm_log.reset_index().rename(columns={"index":"Sample"}))
    labeled.set_index("Sample",inplace=True)
    # drop rows containing NaN values
    labeled = labeled.dropna(axis=0)
    # return dataframe
    return labeled

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
    if balance=="Up":
        sm = SMOTE(random_state=42)
        train_X, y_train = sm.fit_resample(train_X,y_train)
    # for X_train and X_test: scale data to a z-score
    scalar = StandardScaler()
    X_train = scalar.fit_transform(train_X)
    X_test = scalar.fit_transform(test_X)
    # return training and test data
    return X_train, y_train, X_test, y_test

def rf_iterative(stressor,tpmdf,fidf,max_important_features,stepsize,hypdict):
    # subset the feature importance dataframe to the top n features, where n=max_important_features (int)
    topmax = fidf.head(n=max_important_features)

    # initiate classifier
    rfc = RandomForestClassifier(bootstrap=hypdict["bootstrap"],
                                ccp_alpha=hypdict["ccp_alpha"],
                                class_weight=hypdict["class_weight"],
                                criterion=hypdict["criterion"],
                                max_depth=hypdict["max_depth"],
                                max_features=hypdict["max_features"],
                                max_leaf_nodes=hypdict["max_leaf_nodes"],
                                max_samples=hypdict["max_samples"],
                                min_impurity_decrease=hypdict["min_impurity_decrease"],
                                min_samples_leaf=hypdict["min_samples_leaf"],
                                min_samples_split=hypdict["min_samples_split"],
                                min_weight_fraction_leaf=hypdict["min_weight_fraction_leaf"],
                                n_estimators=hypdict["n_estimators"],
                                n_jobs=hypdict["n_jobs"],
                                oob_score=hypdict["oob_score"],
                                random_state=hypdict["random_state"],
                                verbose=hypdict["verbose"],
                                warm_start=hypdict["warm_start"])

    # run the random forest iteratively, saving the scores each time
    setsizes = []
    accuracy = []
    f1_class0 = []
    f1_class1 = []
    for i in range(stepsize,max_important_features,stepsize):
        # append set size to list
        setsizes.append(max_important_features-i)
        # pull out the i least important genes
        leastgenes = list(topmax.tail(n=i)["GeneID"])
        # drop these genes from the TPM dataframe (not yet transformed)
        subtpm = tpmdf.drop(leastgenes,axis=1)
        # transform the data
        transtpm = pre_split_transform(subtpm,False,False)
        # split the train and test sets
        X_train, y_train, X_test, y_test = split_prep_stressor(stressor,transtpm,"Up")
        # fit model on training data
        rfc.fit(X_train,y_train)
        # make predictions on test set
        y_pred = rfc.predict(X_test)
        # generate scores and save into lists
        accuracy.append(accuracy_score(y_test,y_pred))
        f1 = list(f1_score(y_test,y_pred,average=None))
        f1_class0.append(f1[0])
        f1_class1.append(f1[1])

    # stick lists together into a dataframe
    setsize_scores = pd.DataFrame(list(zip(setsizes,accuracy,f1_class0,f1_class1)),columns=["Set Size","Accuracy","F1 Control","F1 Stressed"])
    return setsize_scores
