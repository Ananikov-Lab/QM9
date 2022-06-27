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


def generate_embeddings_dict(embeds_reag, embeds_prod, smiles_reag, smiles_prod):
    embeddings_of_smiles = {}

    all_smiles_reag = []
    all_smiles_prod = []

    embeddings_reag = np.load(embeds_reag)
    embeddings_prod = np.load(embeds_prod)

    with open(smiles_reag, 'rt') as f:
        for cmpd in f:
            all_smiles_reag.append(cmpd.rstrip('\n'))

    with open(smiles_prod, 'rt') as f:
        for cmpd in f:
            all_smiles_prod.append(cmpd.rstrip('\n'))

    for i in range(len(all_smiles_reag)):
        embeddings_of_smiles[all_smiles_reag[i]] = embeddings_reag[i]
    for i in range(len(all_smiles_prod)):
        embeddings_of_smiles[all_smiles_prod[i]] = embeddings_prod[i]

    return embeddings_of_smiles


def reduct(method_name, n_components, embeddings, random_state=42):
    function = globals()[method_name]
    X = embeddings
    X_r = function(n_components, random_state).fit_transform(X)
    
    return X_r


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('embeds_reag_path', type=str)
    parser.add_argument('embeds_prod_path', type=str)
    parser.add_argument('smiles_reag_path', type=str)
    parser.add_argument('smiles_prod_path', type=str)
    parser.add_argument('reactions_path', type=str)
    parser.add_argument('method_name', type=str)
    parser.add_argument('n_components', type=str)
    parser.add_argument('random_state', type=str)
    parser.add_argument('output_file', type=str)
    args = parser.parse_args()

    embeddings_of_reacts = []

    embeddings_of_smiles = generate_embeddings_dict(
        args.embeds_reag_path, args.embeds_prod_path, args.smiles_reag_path, args.smiles_prod_path)
    with open(args.reactions_path, 'rb') as f:
        reactions = pkl.load(f)

    for reacts in reactions.mpairs_of_reacts:
        embedding_react = embeddings_of_smiles[Chem.MolToSmiles(reacts.reag[0].mol_structure)] \
                          - embeddings_of_smiles[Chem.MolToSmiles(reacts.prod[0].mol_structure)]
        embeddings_of_reacts.append(embedding_react)

    X_r = reduct(args.method_name, args.n_components, args.embeddings_of_reacts, args.random_state)

    with open(args.output_file, 'wb') as f:
        pkl.dump(X_r, f)
