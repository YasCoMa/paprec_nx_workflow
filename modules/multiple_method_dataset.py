#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 22:39:05 2022

@author: yasmmin
"""

import os
import re
import urllib.request
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
h = .02  # step size in the mesh

from sklearn import svm
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import AdaBoostClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.naive_bayes import BernoulliNB
from sklearn import tree
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix, roc_auc_score, auc, precision_recall_curve
from sklearn.metrics import make_scorer

from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier

import statistics as st

def scoring_prauc(ytrue, ypred, **kwargs):
    prec, rec, _ = precision_recall_curve(ytrue, ypred)
    prauc = auc(prec, rec)
    return prauc

class Implementation_vaxijen:
    descriptors={}
    
    def __init__(self):
        self.descriptors=self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open("e_table_descriptors.tsv","r")
        for line in f:
            l=line.replace("\n","").split(" ")
            aa=l[0]
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
    
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l.replace('\t','_').replace('.','_')
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')
        
        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )
        
        g.close()
        f.close()
    
    def _calculate_features(self, seq, mode):
        aac={}
        
        n=len(seq)
        
        l=8
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for i in range(5):
                for j in range(5):
                    if(mode=='auto'):
                        cond=(i==j)
                    if(mode=='cross'):
                        cond=(i!=j)
                        
                    if(cond):
                        s=0
                        
                        k=0
                        while k<(n-l):
                            aa=seq[k]
                            flag=False
                            if(aa in self.descriptors.keys()):
                                e=self.descriptors[aa]
                                flag=True
                            
                            aal=seq[k+l]
                            flag2=False
                            if(aal in self.descriptors.keys()):
                                el=self.descriptors[aal]
                                flag2=True
                            
                            if(flag and flag2):
                                s+= (e[j]*el[j]) / (n-l)
                            
                            k+=1
                        
                        aac["l"+str(l)].append(s)
        return aac                
                        
        
    def build_dataset_matrix_variance(self, mode, ide, method):
        fasta_pos=ide+"/"+method+"/dataset_pos.fasta"
        fasta_neg=ide+"/"+method+"/dataset_neg.fasta"
        
        classes={"pos": 1, "neg": 0}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 26):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, mode)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()
        
class Implementation_vaxijenModified:
    descriptors={}
    
    def __init__(self):
        self.descriptors=self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open("descriptors_pmc5549711.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            aa=l[0]
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
    
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')
        
        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )
        
        g.close()
        f.close()
    
    def build_fastas(self):
        i=0
        g=open("alternative_dataset/dataset_pos.fasta","w")
        f=open("alternative_dataset/positive_set.tsv","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            if(i>0):
                g.write(">%s\n%s\n" %( l[0]+"_"+l[2], l[1] ) )
            i+=1
        f.close()
        g.close()
        
        i=0
        g=open("alternative_dataset/dataset_neg.fasta","w")
        f=open("alternative_dataset/negative_set.tsv","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            if(i>0):
                g.write(">%s\n%s\n" %( l[0]+"_"+l[2], l[1] ) )
            i+=1
        f.close()
        g.close()
        
    def _get_optimal_lag(self, fasta_pos, fasta_neg):
        seqs=set()
        f=open(fasta_pos, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        f=open(fasta_neg, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        n=len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1    
        if(i<4):
            i=4
        return i
    
    def _calculate_features(self, seq, max_lag):
        aac={}
        
        n=len(seq)
        
        l=max_lag
        
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for j in range(6):
                mean_=0
                for aa in seq:
                    el=self.descriptors[aa]
                    mean_+=el[j]
                mean_=mean_/n
                
                s=0
                i=0
                while i<(n-l):
                    aa=seq[i]
                    flag=False
                    if(aa in self.descriptors.keys()):
                        e=self.descriptors[aa]
                        va=e[j] - mean_
                        flag=True
                    
                    aal=seq[i+l]
                    flag2=False
                    if(aal in self.descriptors.keys()):
                        el=self.descriptors[aal]
                        vb=el[j] - mean_
                        flag2=True
                    
                    if(flag and flag2):
                        s+= (va * vb) / (n-l)
                    
                    i+=1
                        
                aac["l"+str(l)].append(s)
        return aac                
                        
        
    def build_dataset_matrix_variance(self, mode, ide, method):
        fasta_pos=ide+"/"+method+"/dataset_pos.fasta"
        fasta_neg=ide+"/"+method+"/dataset_neg.fasta"
        
        max_lag=self._get_optimal_lag(fasta_pos, fasta_neg)
        print("max lag", max_lag)
        
        classes={"pos": 1, "neg": 0}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 7):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, max_lag)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace(' ','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace(' ','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()
 
class Implementation_new_158:
    descriptors={}
    
    def __init__(self):
        self.descriptors=self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        df=pd.read_csv('filtered_aaindex.tsv', sep='\t')
        for k in df.columns:
            if(k!='id' and k!='name' and k!='description'):
                descriptors[k]=list(df[k])
        
        return descriptors
    
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l.replace('\t','_').replace('.','_')
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')
        
        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )
        
        g.close()
        f.close()
    
    def _calculate_features(self, seq, mode):
        aac={}
        
        n=len(seq)
        
        l=1
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for i in range(len(self.descriptors['M'])):
                init=0
                for s in seq:
                    init+=self.descriptors[s][i]
                aac["l"+str(l)].append(init)
        return aac                
                        
        
    def build_dataset_matrix_variance(self, mode, ide, method):
        fasta_pos=ide+"/"+method+"/dataset_pos.fasta"
        fasta_neg=ide+"/"+method+"/dataset_neg.fasta"
        
        classes={"pos": 1, "neg": 0}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 26):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, mode)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()
         
class Screening_classifier:
        
    def _get_optimal_lag(self, fasta_pos, fasta_neg):
        seqs=set()
        f=open(fasta_pos, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        f=open(fasta_neg, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        n=len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1        
        return i
    
    def prepare_random_datasets(self, ide, met):
        modes =['auto','cross']
        maxi_lag=8
        if(met=='method2' or met=='method3'):
            modes =['auto']
            fasta_pos=ide+"/"+met+"/dataset_pos.fasta"
            fasta_neg=ide+"/"+met+"/dataset_neg.fasta"
            maxi_lag=self._get_optimal_lag(fasta_pos, fasta_neg)
            
        for i in range(5):
            for m in modes:
                df=pd.read_csv(ide+"/"+met+"/"+m+"_dataset.tsv", sep="\t")
                pieces=[]
                for lag in range(1, maxi_lag+1):
                    auxdf=df[ (df['class']==1) & (df['lag']=='l'+str(lag)) ]
                    if(len(auxdf)>0):
                        ns=100
                        if(len(auxdf)<ns):
                            ns=len(auxdf)
                        chosen=random.sample( list(auxdf.index), ns)
                        pieces.append(df.iloc[chosen, :])
                    
                    auxdf=df[ (df['class']==0) & (df['lag']=='l'+str(lag)) ]
                    if(len(auxdf)>0):
                        ns=100
                        if(len(auxdf)<ns):
                            ns=len(auxdf)
                        chosen=random.sample( list(auxdf.index), ns)
                        pieces.append(df.iloc[chosen, :])
                        
                if(len(pieces)>0):    
                    dfv=pd.concat( pieces )
                    dfv.to_csv(ide+"/"+met+"/"+"random_datasets/"+str(i)+"_"+m+"_dataset.tsv", sep="\t", index=False)
                
    def test_cross_validation(self, ide, met):
        kernel = 1.0 * RBF(1.0)
        
        clfs={'svm-rbf': "svm.SVC(kernel='rbf')", 
              'sgd': 'SGDClassifier(loss="hinge", penalty="l2", max_iter=5)',
              'nearestCentroid': 'NearestCentroid()',
              'adaboost': 'AdaBoostClassifier()',
              'gaussianProcess': 'GaussianProcessClassifier(kernel=kernel, random_state=0)',
              'BernoulliNaiveBayes': 'BernoulliNB()',
              'DecisionTree': 'tree.DecisionTreeClassifier()',
              'GradientBoosting': "GradientBoostingClassifier( n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)",
              'NeuralNetwork': "MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)",
              'randomForest': 'RandomForestClassifier(max_depth=5)'
        }
        
        f=open(ide+"/"+met+"/"+"result_cross-validation.txt","w")
        f.write("dataset;mode;lag;classifier;f1;precision;recall;accuracy;roc_auc;pr_auc\n")
        modes =['auto','cross']
        maxi_lag=8
        if(met=='method2' or met=='method3'):
            modes =['auto']
            fasta_pos=ide+"/"+met+"/dataset_pos.fasta"
            fasta_neg=ide+"/"+met+"/dataset_neg.fasta"
            maxi_lag=self._get_optimal_lag(fasta_pos, fasta_neg)
        
        prauc_scorer = make_scorer(scoring_prauc, greater_is_better=True)
            
        for ds in range(5):
            for m in modes:
                for lag in range(1, maxi_lag+1):
                    
                    df=pd.read_csv(ide+"/"+met+"/"+"random_datasets/"+str(ds)+"_"+m+"_dataset.tsv", sep="\t")
                    x = df[ df['lag']=='l'+str(lag) ]
                    if(len(x)>1):
                        y=x.iloc[:,1]
                        X=[]
                        for i in range(len(x)):
                            aux=[]
                            values=x.iloc[i,3].split(",")
                            for v in values:
                                aux.append(float(v))
                            X.append(aux)
                            
                        for classifier in clfs:
                            clf = eval(clfs[classifier])
                            
                            f1=cross_val_score(clf, X, y, scoring='f1', cv=10)
                            precision=cross_val_score(clf, X, y, scoring='precision', cv=10)
                            recall=cross_val_score(clf, X, y, scoring='recall', cv=10) 
                            accuracy=cross_val_score(clf, X, y, scoring='accuracy', cv=10)
                            roc_aucs=cross_val_score(clf, X, y, scoring='roc_auc', cv=10)
                            praucs=cross_val_score(clf, X, y, scoring=prauc_scorer, cv=10)
                            
                            for i in range(len(f1)):
                                f.write(str(ds)+";"+m+";l"+str(lag)+";"+classifier+";"+str(f1[i])+";"+str(precision[i])+";"+str(recall[i])+";"+str(accuracy[i])+";"+str(roc_aucs[i])+";"+str(praucs[i])+"\n")
                            
                            ides="ds"+str(ds)+"_"+m+"_l"+str(lag)+"_"+classifier
                            self._auc_roc(ides, X, y, clf, ide, met)
        f.close()
        
    def _auc_roc(self, ides, X, y, clf, ide, met):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)
        clf.fit(X_train, y_train)
        try:
            predictions = clf.predict_proba(X_test)
            predictions = np.array(predictions)[:, 1]
        except:
            predictions = clf.predict(X_test)
            
        false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, predictions)
        roc_auc = auc(false_positive_rate, true_positive_rate)
        if(roc_auc >= 0.6):
            f = plt.figure()
            plt.title('ROC Curve - Compilation)')
            plt.plot(false_positive_rate, true_positive_rate, 'b', label='AUC = %0.2f'% roc_auc)
            plt.legend(loc='lower right')
            plt.plot([0,1],[0,1],'r--')
            plt.xlim([-0.1,1.2])
            plt.ylim([-0.1,1.2])
            plt.ylabel('True Positive Rate')
            plt.xlabel('False Positive Rate')
            f.savefig(ide+"/"+met+"/"+'performance_aucroc/'+ides+'roc_auc.png')        
        
    def check_standard_deviation(self, ide, met):
        g=open(ide+"/"+met+"/"+"summary_cross_validation_result.tsv","w")
        g.write("dataset\tmode\tlag\tclassifier\tmean f1\tst dev f1\tmean precision\tst dev precision\tmean recall\tst dev recall\tmean accuracy\tst dev accuracy\tmean roc_auc\tst dev roc_auc\tmean pr_auc\tst dev pr_auc\n")
        #g.write("dataset\tmode\tlag\tclassifier\tmean f1\tst dev f1\tmean precision\tst dev precision\tmean recall\tst dev recall\tmean accuracy\tst dev accuracy\tmean roc_auc\tst dev roc_auc\n")
        c=0
        
        ant=""
        f1=[]
        acc=[]
        rec=[]
        prec=[]
        rocauc=[]
        prauc=[]
        f=open(ide+"/"+met+"/"+"result_cross-validation.txt","r")
        for line in f:
            l=line.replace("\n","").split(";")
            if(c>0):
                if('\t'.join(l[:4]) != ant):
                    if(ant!=""):
                        #g.write(ant+"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" %(st.mean(f1), st.stdev(f1), st.mean(prec), st.stdev(prec), st.mean(rec), st.stdev(rec), st.mean(acc), st.stdev(acc), st.mean(rocauc), st.stdev(rocauc) ) )
                        g.write(ant+"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" %(st.mean(f1), st.stdev(f1), st.mean(prec), st.stdev(prec), st.mean(rec), st.stdev(rec), st.mean(acc), st.stdev(acc), st.mean(rocauc), st.stdev(rocauc), st.mean(prauc), st.stdev(prauc) ) )
                    ant='\t'.join(l[:4])
                    f1=[]
                    acc=[]
                    rec=[]
                    prec=[]
                    rocauc=[]
                    prauc=[]
                
                f1.append(float(l[4]))
                prec.append(float(l[5]))
                rec.append(float(l[6]))
                acc.append(float(l[7]))
                rocauc.append(float(l[8]))
                prauc.append( float(l[9]) )
            c+=1
        f.close()
        
        if(ant!=""):
            #g.write(ant+"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" %(st.mean(f1), st.stdev(f1), st.mean(prec), st.stdev(prec), st.mean(rec), st.stdev(rec), st.mean(acc), st.stdev(acc), st.mean(rocauc), st.stdev(rocauc) ) )
            g.write(ant+"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" %(st.mean(f1), st.stdev(f1), st.mean(prec), st.stdev(prec), st.mean(rec), st.stdev(rec), st.mean(acc), st.stdev(acc), st.mean(rocauc), st.stdev(rocauc), st.mean(prauc), st.stdev(prauc) ) )
        g.close()

class Test_feature_selection:
    
    def prepare_condensed_dataset(self, ides, met):
        
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        for m in modes:
            flag_stop=False
            f=open(ides+"/"+met+"/"+m+"_mapping_name_features.tsv", "w")
            f.close()
        
            index=1
            
            df=pd.read_csv(ides+"/"+met+"/"+m+"_dataset.tsv", sep="\t")
            y=[]
            auxx={}
            aux=[]
            for i in range(len(df)):
                ide=df.iloc[i, 0]+'-'+str(df.iloc[i, 1])
                ide=ide.replace('\t','').replace(' ','').replace(':','')
                lag=df.iloc[i,2]
                
                if(not ide in auxx.keys()):
                    if(len(auxx.keys())>0):
                        flag_stop=True
                        
                    auxx[ide]=[]
                    y.append(df.iloc[i, 1])
                
                h=1
                aux=[]
                values=df.iloc[i,3].split(",")
                for v in values:
                    aux.append(str(v))
                    
                    if(not flag_stop):
                        with open(ides+"/"+met+"/"+m+"_mapping_name_features.tsv", "a") as f:
                            f.write("%s\t%s\n" %("f"+str(index), lag+"_f"+str(h)) )
                    h+=1
                    index+=1
                
                auxx[ide]+=aux
                
            f=open(ides+"/"+met+"/"+m+"_dataset_consolidated.tsv", "w")
            f.write("protein\tclass\tfeatures\n")
            for k in auxx.keys():
                info=k.split("-")
                f.write( ('\t'.join(info))+"\t"+(','.join(auxx[k]))+"\n")
            f.close()
            
    def prepare_random_datasets(self, ide, met):
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        for i in range(5):
            for m in modes:
                #print(ide+"/"+met+"/"+m+"_dataset_consolidated.tsv")
                df=pd.read_csv(ide+"/"+met+"/"+m+"_dataset_consolidated.tsv", sep="\t")
                pieces=[]
                auxdf=df[ (df['class']==1) ]
                ns=100
                if(len(auxdf)<ns):
                    ns=len(auxdf)
                chosen=random.sample( list(auxdf.index), ns)
                pieces.append(df.iloc[chosen, :])
                
                auxdf=df[ (df['class']==0) ]
                ns=100
                if(len(auxdf)<ns):
                    ns=len(auxdf)
                chosen=random.sample( list(auxdf.index), ns)
                pieces.append(df.iloc[chosen, :])
                    
                dfv=pd.concat( pieces )
                dfv.to_csv(ide+"/"+met+"/"+"random_datasets/"+str(i)+"_"+m+"_dataset_condensed.tsv", sep="\t", index=False)
                
    def execute_feature_selection(self, ide, met):
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        
        mapp={}
        for m in modes:
            mapp[m]={}
            f=open(ide+"/"+met+"/"+m+"_mapping_name_features.tsv", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                mapp[m][l[0]]=l[1]
            f.close()
        
        clfs={'randomForest': 'RandomForestClassifier(max_depth=5)',
              'GradientBoosting': "GradientBoostingClassifier( n_estimators=100, learning_rate=1.0, max_depth=5, random_state=0)"
        }
        
        X=[]
        Y=[]
        fg=open(ide+"/"+met+"/"+"feature_selection_results/model_results.tsv","w")
        fg.write("dataset\tselected_features\tmode\tclassifier\taccuracy\tf1\tprecision\trecall\troc_auc\tpr_auc\n")
        
        for ds in range(5):
            for m in modes:
                df=pd.read_csv(ide+"/"+met+"/"+"random_datasets/"+str(ds)+"_"+m+"_dataset_condensed.tsv", sep="\t")
                Y=df.iloc[:,1]
                X=[]
                feas=[]
                for i in range(len(df)):
                    aux=[]
                    values=df.iloc[i,2].split(",")
                    j=1
                    for v in values:
                        aux.append(float(v))
                        if(i==0):
                            feas.append("f"+str(j))
                        j+=1
                    X.append(aux)
                
                X=pd.DataFrame(X)
                X.columns=feas
                X.to_csv(ide+"/"+met+"/"+"feature_selection_results/"+str(ds)+"_"+m+"_features_selection_dataset.tsv", sep='\t')
                
                for classifier in clfs.keys():
                    #print("-----", ds, m, classifier)
                    
                    #forest = RandomForestClassifier(n_jobs=-1,  max_depth=5)
                    forest=eval(clfs[classifier])
                    forest.fit(X, Y)
                    feat_selector = BorutaPy(forest, n_estimators='auto', random_state=1)
                    feat_selector.fit(X.to_numpy(), Y)
                    
                    feature_ranks = list(zip(X.columns, feat_selector.ranking_, feat_selector.support_))
                    sels=[]
                    f=open(ide+"/"+met+"/"+"feature_selection_results/"+classifier+"_"+str(ds)+"_"+m+"_features_selection_result.tsv","w")
                    # iterate through and print out the results
                    for feat in feature_ranks:
                        f.write('Feature: {:<25} Rank: {},  Keep: {} \n'.format(feat[0], feat[1], feat[2]))
                        if(feat[2]):
                            sels.append(mapp[m][feat[0]])
                    f.close()
                        
                    X_filtered = feat_selector.transform(X.to_numpy())
                    if( len(X_filtered[0,:]) > 0):
                        X_test=X_filtered[:20000,:]
                        y_test=Y[:20000]
                        clf = AdaBoostClassifier()
                        clf.fit(X_filtered, Y)
                        predictions = clf.predict(X_test)
                        try:
                            predprob = clf.predict_proba(X_test)
                            predprob = np.array(predprob)[:, 1]
                        except:
                            predprob = predictions
                        
                        roc_auc = roc_auc_score(y_test, predprob)
                        precision, recall, _ = precision_recall_curve(y_test, predprob)
                        pr_auc = auc(recall, precision) 
                        
                        fg.write("%i\t%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %(ds, (','.join(sels)), m, classifier, accuracy_score(y_test, predictions), f1_score(y_test, predictions), precision_score(y_test, predictions), recall_score(y_test, predictions), roc_auc, pr_auc ) )
                        
                        """with open("feature_selection_results/"+classifier+"_"+str(ds)+"_"+m+"_model_result.tsv","w") as fg:
                            fg.write('accuracy:'+str(accuracy_score(y_test, predictions))+"\n")
                            fg.write('f1_score:'+str(f1_score(y_test, predictions))+"\n")
                            fg.write('recall:'+str(recall_score(y_test, predictions))+"\n")
                            fg.write('precision:'+str(precision_score(y_test, predictions))+"\n")"""
        fg.close()
        
        
class Running_config:
    def _init_folders_for_new_ds(self, identificador, m):
        os.system("mkdir "+identificador+"/"+m)
        os.system("cp "+identificador+"/dataset_* "+identificador+"/"+m+"/")
        os.system("mkdir "+identificador+"/"+m+"/random_datasets")
        os.system("mkdir "+identificador+"/"+m+"/performance_aucroc")
        os.system("mkdir "+identificador+"/"+m+"/feature_selection_results")
            
    def run(self):
        met1=Implementation_vaxijen()
        met2=Implementation_vaxijenModified()
        met3=Implementation_new_158()
        
        datasets=['bcipep_dataset','hla_dataset', 'gram+_dataset','gram-_dataset']
        #datasets=[ 'gram+_dataset','gram-_dataset']
        #datasets=[ 'gram-_dataset']
        #datasets=['vipr_db']
        
        methods=['method1', 'method2', 'method3']
        methods=['method3']
        methods=['method1', 'method2']
        
        for ds in datasets:
            for m in methods:
                if(not os.path.isdir(ds+"/"+m)):
                    self._init_folders_for_new_ds(ds, m)
                    
                print("\n>>>>>>>>>>>>>>>>>> Stage ", ds, m, "<<<<<<<<<<<<<<<<<<<<<<<")
                    
                print("----------- Preparing dataset")
                modes =['auto','cross']
                instance=met1
                if(m=='method2'):
                    modes = ['auto']
                    instance=met2
                if(m=='method3'):
                    modes = ['auto']
                    instance=met3
                    
                for mo in modes:
                    instance.build_dataset_matrix_variance(mo, ds, m)  
                    
                print("----------- Screening_classifier")
                b=Screening_classifier()
                print("\texecuting 1...")
                b.prepare_random_datasets(ds, m)
                print("\texecuting 2...")
                b.test_cross_validation(ds, m)
                print("\texecuting 3...")
                b.check_standard_deviation(ds, m)
                    
                print("----------- Test_feature_selection")
                c=Test_feature_selection()
                print("\texecuting 1...")
                c.prepare_condensed_dataset(ds, m)
                print("\texecuting 2...")
                c.prepare_random_datasets(ds, m)
                print("\texecuting 3...")
                c.execute_feature_selection(ds, m)
                
r = Running_config()
r.run()
