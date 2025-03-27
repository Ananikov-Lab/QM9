from argparse import ArgumentParser
import pickle as pkl
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt


def clustering(method_name, X_r, parameters={}):
    function = {
        'KMeans': KMeans,
        'SpectralClustering': SpectralClustering,
        'AgglomerativeClustering': AgglomerativeClustering
    }

    if method_name not in function.keys():
        raise ValueError
    else:
        method = function[method_name]

    sc_model = method().set_params(**parameters).fit(X_r)

    return sc_model


def metric(model, X_r, metric_='euclidean'):
    labels = model.labels_
    metric_value = metrics.silhouette_score(X_r, labels, metric=metric_)
    print('Your metric:', metric_value)

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
    parser.add_argument('--input', dest='input_embeds', type=str,
                        help='file with embeds of reactions after dim_redc')
    parser.add_argument('--method', dest='method_name', type=str,
                        help='name of clustering method')
    parser.add_argument('--metric', dest='metric_name', type=str,
                        help='name of prefer metric')
    parser.add_argument('--plot', dest='path_to_plot', type=str,
                        help='path of plot`s figure')
    parser.add_argument('--model', dest='path_save_model', type=str,
                        help='path of save model')
    parsed, unknown = parser.parse_known_args()
    params = {}

    for arg in unknown:
        if arg.startswith('-'):
            parser.add_argument(arg.split('=')[0], dest=arg.split('=')[0][1:], type=int)
            params[arg.split('=')[0][1:]] = ''

    args = parser.parse_args()

    args_param = vars(args)

    for name in params.keys():
        params[name] = args_param[name]
        
    with open(args.input_embeds, 'rb') as f:
        X_r = pkl.load(f)

    sc_model = clustering(args.method_name, X_r, parameters=params)

    with open(args.path_save_model, 'wb') as f:
        pkl.dump(sc_model, f)

    metric = metric(sc_model, X_r, args.metric_name)

    vizualization(sc_model, X_r, args.path_to_plot)
