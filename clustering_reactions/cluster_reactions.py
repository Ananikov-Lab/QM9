from argparse import ArgumentParser
import pickle as pkl
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn import metrics
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
import seaborn as sns
import umap
import umap.plot
import gzip
import os
import re
import matplotlib.pyplot as plt


def clusterization(method_name, n_clusters, X_r):
    function = globals()[method_name]
    sc_model = function(n_clusters).fit(X_r)

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
    parser.add_argument('input_model', type=str)
    parser.add_argument('method_name', type=str)
    parser.add_argument('n_clusters', type=str)
    parser.add_argument('metric_name', type=str)
    parser.add_argument('path_to_plot', type=str)
    args = parser.parse_args()

    with open(args.input_model, 'rb') as f:
        X_r = pkl.load(f)

    sc_model = clusterization(args.method_name, args.n_clusters, X_r)

    metric = metric(sc_model, X_r, args.metric_name)

    vizualization(sc_model, X_r, args.path_to_plot)
