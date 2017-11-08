#!/usr/bin/env python

import sys
import csv
import os
import time
import psutil

import numpy as np
import subprocess
from collections import Counter
from scipy.stats import fisher_exact

DATA_DIR = 'data'
CROW_HOME = 'crow'

### Utility functions ###

def load_csv(filename, delimiter=',', comment='#'):
    # loads csv file into array of rows (arrays)
    csv.field_size_limit(sys.maxsize)
    fp = open(filename)
    if not fp:
        raise Exception("Error: Cannot open file: %s" % filename)
    
    reader = csv.reader(fp, delimiter=delimiter)
    data = []
    for row in reader:
        if row[0][0] == comment:
            continue
        data.append(row)
    fp.close()
    return data

def load_numpy(filename):
    if os.path.isfile(filename) == False:
        print "Error: Cannot open file: %s" % filename
        return None
    d = np.load(filename)
    if 'indices' in d:
        # sparse
        X = csr_matrix((d['data'], d['indices'], d['indptr']), shape=d['shape'])
        return X
    elif 'data' in d:
        return d['data']
    else:
        print "Error: No numpy data in file %s" % filename
        return None

### Path converter ###

def file_concat(parts):
    path = ''
    for i in range(len(parts)):
        path = os.path.join(path, parts[i])
    return path

def to_path(*parts):
    return file_concat(parts)
    

### Analyzer functions ###

def vectorize(X):
    # Determine cluster ID according to max value
    R = np.zeros((X.shape[0],))
    for i in range(X.shape[0]):
        _max = 0
        jmax = 0
        for j in range(X.shape[1]):
            if X[i,j] > _max:
                jmax = j
                _max = X[i,j]
        R[i] = jmax
    return R

def enrichment_analysis(class_items):
    # Calculate fisher exact test of each class of each cluster
    # Contingency matrix 2x2 with the following frequencies:
    # [ ['This cluster & this class', 'other classes & this cluster' ], 
    #   ['other clusters & this class', 'other clusters & other classes' ] ]
    clusters = {}
    for key, value in class_items:
        for i in value:
            if i not in clusters:
                clusters[i] = []
            clusters[i].append(key)
    
    sum_class = {}
    sum_cluster = {}
    for value in sorted(clusters.keys()):
        d = dict(Counter(clusters[value]))
        sum_cluster[value] = sum(d.values())
        for v in d:
            if v not in sum_class:
                sum_class[v] = 0
            sum_class[v] += d[v]
    
    enrichments = {}
    for cluster in sorted(clusters.keys()):
        d = dict(Counter(clusters[cluster]))
        cluster_enrichments = {}
        sum_other = sum([sum_cluster[v] for v in clusters.keys() if v != cluster])
        for cls in d:
            oc = sum_class[cls] - d[cls]
            oclusters = sum_cluster[cluster] - d[cls]
            C = np.array([[d[cls], oclusters], [oc, sum_other - oc]])
            cluster_enrichments[cls] = fisher_exact(C)

        enrichments[cluster] = cluster_enrichments
    
    return enrichments

### Helper functions ###

def sort_by_pvalue(d):
    return sorted(d.items(), key=lambda x: x[1][1])

def unique(lst):
    return sorted(list(set(lst)))

def eformat(f, prec, exp_digits):
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return "%s \\cdot 10^{ %+0*d }"%(mantissa, exp_digits+1, int(exp))

### Print functions ###

def print_interactions(S, row_clusters, column_clusters, limit=10):
    slist = []
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            slist.append((i, j, S[i,j]))
    
    sorted_slist = sorted(slist, key=lambda x: x[2], reverse=True)
    for x in range(limit):
        i, j, val = sorted_slist[x]
        print "Interaction %d" % x, val
        a = row_clusters[i]
        b = column_clusters[j]
        
        print "Cancer cluster %02d" % i, [x[0] for x in sort_by_pvalue(a)[:3]]
        print "GO cluster     %02d" % j, [x[0] for x in sort_by_pvalue(b)[:3]]
        print


def visu_enrichment(enrichment_dict, use_latex=False):
    for value in enrichment_dict:
        d = enrichment_dict[value]
        
        print "Cluster %d:" % value
        for x in sort_by_pvalue(d)[:5]:
            enrichment = '%.2f' % x[1][0]
            disease = x[0]
            print '%-6s' % enrichment, ' p-value: %-20s' % str(x[1][1]), disease
        print
    
    return enrichment_dict

### Call factorization ###

def factorize(k1=10, k2=10):
    data_src = to_path(DATA_DIR, 'TCGA-Methyl-cancer.npz')
    data_dst = to_path(CROW_HOME, 'data', 'TCGA-Methyl-cancer.npz')
    
    if not os.path.exists(data_dst):
        shutil.copy(data_src, data_dst)
    
    n_cpu = psutil.cpu_count()
    blocks = '1x%d' % n_cpu
    cmd = 'crow -b %s -i 1000 -k1 %d -k2 %d TCGA-Methyl-cancer.npz' % (blocks, k1, k2)
    cmd = cmd.split(' ')
    print "Factorization started, please wait..."
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    print out
    if err != '':
        print err
        raise Exception("Error encountered during factorization")

### Entry point ###

def main():
    rows = load_csv(to_path(DATA_DIR, 'row-labels.csv'))
    cols = load_csv(to_path(DATA_DIR, 'col-labels.csv'))
    
    ftype = 'long'
    k1 = 25
    k2 = 30
    
    factorize(k1=25, k2=30)
    
    dirname = to_path(CROW_HOME, 'results')
    if not os.path.isdir(dirname):
        raise Exception("Directory not found: %s" % dirname)
    U = load_numpy('%s/U.npz' % dirname)
    S = load_numpy('%s/S.npz' % dirname)
    V = load_numpy('%s/V.npz' % dirname)

    U = vectorize(U)
    V = vectorize(V)

    cancer_classes = {}
    for i, (row, cancer) in enumerate(rows):
        if cancer not in cancer_classes:
            cancer_classes[cancer] = []
        cancer_classes[cancer].append(U[i])
    
    gene_classes = {}
    for i, (col, gene) in enumerate(cols):
        if gene != '':
            for name in gene.split(';'):
                if name not in gene_classes:
                    gene_classes[name] = []
                gene_classes[name].append(V[i])

    print "---------------------------------------"
    print "Row clusters"
    row_clusters = visu_enrichment(enrichment_analysis(cancer_classes.items()))
    print "---------------------------------------"
    print "Column clusters"
    column_clusters = visu_enrichment(enrichment_analysis(gene_classes.items()))
    print "---------------------------------------"
    print "10 Strongest Interactions in S matrix:"
    print_interactions(S, row_clusters, column_clusters, limit=10)

    
    
if __name__ == '__main__':
    if 'CROW_HOME' in os.environ:
        CROW_HOME = os.environ['CROW_HOME']
    dcomp_file = to_path(CROW_HOME, 'docker-compose.yml')
    if not os.path.exists(dcomp_file):
        print "You may need to install CROW framework."
        print "%s missing. Is CROW_HOME set?" % dcomp_file
    else:
        main()