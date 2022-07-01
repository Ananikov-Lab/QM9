from argparse import ArgumentParser
import pickle as pkl
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn import metrics
import pandas as pd
import numpy as np
import seaborn as sns
import gzip
import os
import re
import matplotlib.pyplot as plt


def clustering(method_name, n_clusters, X_r, parameters={}):
    function = {
                'KMeans': KMeans,
                'SpectralClustering': SpectralCluastering,
                'AgglomerativeClustering': AgglomerativeClustering
    }
    
    if method_name not in function.keys():
        raise ValueError
    else:
        method = function[method_name]
    
    if len(parameters.keys() - method.__dict__.keys()) > 0:
        raise ImportError
    
    
    sc_model = method(n_clusters, parameters).fit(X_r)

    return sc_model


def metric(model, X_r, metric_='euclidean'):
    labels = model.labels_
    metric_value = metrics.silhouette_score(X_r, labels, metric_)
    print(metric_value)

    return metric_value


def vizualization(model, X_r, path_to_plot):
    labels = model.labels_

    c = np.array(labels)
    x = X_r[:, 0]
    y = X_r[:, 1]
    plt.figure(figsize=(10, 10))

    sc = plt.scatter(x, y, c=c, s=1)
    plt.colorbar(sc)
    plt.savefig(path_to_plot)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', 'input_embeds', type=str, 
                       help='file with embeds of reactions')
    parser.add_argument('--method', 'method_name', type=str, 
                       help='name of clustering method')
    parser.add_argument('--nclusters', 'n_clusters', type=str, 
                       help='number clusters')
    parser.add_argument('--metric', 'metric_name', type=str, 
                       help='name of prefer metric')
    parser.add_argument('--plot', 'path_to_plot', type=str, 
                       help='path of plot`s figure')
    
    params = {}
    
    for arg in unknown:
        if arg.startswith('--'):
            parser.add_argument(arg.split('=')[0], dest=arg.split("=")[0][2:])
            params[arg.split("=")[0][2:]] = ''
    
    args = parser.parse_args()
    
    args_param = vars(args)
    
    for name in params.keys():
        params[name] = args_param[name]

    with open(args.input_embeds, 'rb') as f:
        X_r = pkl.load(f)

    sc_model = clustering(args.method_name, args.n_clusters, X_r, parameters=params)

    metric = metric(sc_model, X_r, args.metric_name)

    vizualization(sc_model, X_r, args.path_to_plot)
