import pandas as pd
from argparse import ArgumentParser


def extract_energy(path_to_cmpds, dict_cmpds, dict_cmpds_smiles):
    with open(path_to_cmpds) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            index = line[0]
            smiles = line[1]
            gibbs_energy = line[2]
            dict_cmpds[index] = float(gibbs_energy)
            dict_cmpds_smiles[index] = smiles
    return dict_cmpds, dict_cmpds_smiles


def data_to_csv(path_to_reags, path_to_alks, path_to_prods, output_file):
    dict_reags_smiles = {}
    dict_alks_smiles = {}
    dict_prods_smiles = {}

    dict_reags = {}
    dict_alks = {}
    dict_prods = {}

    smiles_reags = []
    smiles_alks = []
    smiles_prods = []
    gibbs_energy_react = []

    dict_reags, dict_reags_smiles = extract_energy(path_to_reags, dict_reags, dict_reags_smiles)
    dict_alks, dict_alks_smiles = extract_energy(path_to_alks, dict_alks, dict_alks_smiles)
    dict_prods, dict_prods_smiles = extract_energy(path_to_prods, dict_prods, dict_prods_smiles)

    for index, energy in dict_prods.items():
        index_list = index.split('_')
        reag = index_list[0]
        alk = index_list[1]
        rxn_gibbs_energy = energy - dict_reags[reag] - dict_alks[alk]
        smiles_reags.append(dict_reags_smiles[reag])
        smiles_alks.append(dict_alks_smiles[alk])
        smiles_prods.append(dict_prods_smiles[index])
        gibbs_energy_react.append(rxn_gibbs_energy)

    df = pd.DataFrame({
        'reags': smiles_reags,
        'alks': smiles_alks,
        'prods': smiles_prods,
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
