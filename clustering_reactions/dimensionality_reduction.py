import pickle as pkl
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
import umap
import umap.plot
from mining_pubchem.find_reactions import *
from mining_pubchem.processing_dataset import *

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


def reduct(method_name, embeddings, parameters={}):
    function = {
        'PCA': PCA,
        't-SNE': TSNE,
        'UMAP': umap.UMAP
    }

    if method_name not in function.keys():
        raise ValueError
    else:
        method = function[method_name]

    if len(parameters.keys() - method.__dict__.keys()) > 0:
        raise ImportError

    X = embeddings

    X_r = method().set_params(**parameters).fit_transform(X)

    return X_r


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--emb-r', dest='embeds_reag_path', type=str,
                        help='binary file of reagents embeddings')
    parser.add_argument('--emb-p', dest='embeds_prod_path', type=str,
                        help='binary file of products embeddings')
    parser.add_argument('--smi-r', dest='smiles_reag_path', type=str,
                        help='file of reagents smiles')
    parser.add_argument('--smi-p', dest='smiles_prod_path', type=str,
                        help='file of products smiles')
    parser.add_argument('--reacts', dest='reactions_path', type=str,
                        help='binary file of class "Reactions" with generated reactions')
    parser.add_argument('--method', dest='method_name', type=str,
                        help='name of dimensionality reduction method')
    parser.add_argument('--output', dest='output_file', type=str,
                        help='path of embeddings after dimensionality reductions')
    parsed, unknown = parser.parse_known_args()

    params = {}

    for arg in unknown:
        if arg.startswith('-'):
            parser.add_argument(arg.split('=')[0], dest=arg.split("=")[0][1:], type=int)
            params[arg.split("=")[0][1:]] = ''

    args = parser.parse_args()

    args_param = vars(args)

    for name in params.keys():
        params[name] = args_param[name]

    embeddings_of_reacts = []

    embeddings_of_smiles = generate_embeddings_dict(
        args.embeds_reag_path, args.embeds_prod_path, args.smiles_reag_path, args.smiles_prod_path)
    with open(args.reactions_path, 'rb') as f:
        reactions = pkl.load(f)

    for reacts in reactions.mpairs_of_reacts:
        embedding_react = embeddings_of_smiles[Chem.MolToSmiles(reacts.reag[0].mol_structure)] \
                          - embeddings_of_smiles[Chem.MolToSmiles(reacts.prod[0].mol_structure)]
        embeddings_of_reacts.append(embedding_react)

    X_r = reduct(args.method_name, embeddings_of_reacts, parameters=params)

    with open(args.output_file, 'wb') as f:
        pkl.dump(X_r, f)
