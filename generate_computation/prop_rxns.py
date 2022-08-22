import pandas as pd
from argparse import ArgumentParser


def data_to_csv(path_to_reags, path_to_alks, path_to_prods, output_file):
    dict_reag_smiles = {}
    dict_alkyne_smiles = {}
    dict_prod_smiles = {}

    dict_reag = {}
    dict_alkyne = {}
    dict_prod = {}

    smiles_reag = []
    smiles_alk = []
    smiles_prod = []
    gibbs_energy_react = []

    with open(path_to_reags) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            index = line[0]
            smiles = line[1]
            gibbs_energy = line[2]
            dict_reag[index] = float(gibbs_energy)
            dict_reag_smiles[index] = smiles

    with open(path_to_alks) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            index = line[0]
            smiles = line[1]
            gibbs_energy = line[2]
            dict_alkyne[index] = float(gibbs_energy)
            dict_alkyne_smiles[index] = smiles

    with open(path_to_prods) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            index = line[0]
            smiles = line[1]
            gibbs_energy = line[2]
            dict_prod[index] = float(gibbs_energy)
            dict_prod_smiles[index] = smiles

    for index, energy in dict_prod.items():
        index_ = index.split('_')
        reag = index_[0]
        alk = index_[1]
        rxn_gibbs_energy = energy - dict_reag[reag] - dict_alkyne[alk]
        smiles_reag.append(dict_reag_smiles[reag])
        smiles_alk.append(dict_alkyne_smiles[alk])
        smiles_prod.append(dict_prod_smiles[index])
        gibbs_energy_react.append(rxn_gibbs_energy)

    df = pd.DataFrame({
        'reag': smiles_reag,
        'alkyne': smiles_alk,
        'prod': smiles_prod,
        'gibbs energy': gibbs_energy_react}
    )

    df.to_csv(output_file)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('path_to_reags', type=str,
                        help="First - path of .txt file with reag`s data")
    parser.add_argument('path_to_alks', type=str,
                        help="Second - path of .txt file with alk`s data")
    parser.add_argument('path_to_prods', type=str,
                        help="Third - path of .txt file with prod`s data")
    parser.add_argument('output_file', type=str,
                        help="Fourth - path of output .csv file")
    args = parser.parse_args()

    data_to_csv(args.path_to_reags, args.path_to_alks,
                args.path_to_prods, args.output_file)
