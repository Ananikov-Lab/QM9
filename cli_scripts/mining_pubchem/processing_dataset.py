import tarfile
from tqdm import tqdm
import matplotlib.pyplot as plt
import gzip
import os
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger
from argparse import ArgumentParser
import pickle as pkl
from rdkit.Chem.rdchem import KekulizeException
import numpy as np

RDLogger.DisableLog('rdApp.*')


class Compound:
    def __init__(self, formula, gibbs_energy, smiles, dataset_presence):
        try:
            self.formula = formula
            self.gibbs_energy = gibbs_energy
            self.smiles = smiles
            self.dataset_presence = dataset_presence
            self.mol_structure = Chem.MolFromSmiles(smiles[0])
            self.inchi = Chem.MolToInchi(self.mol_structure)
        except BaseException:
            self.mol_structure = None


class KiloMolecules:
    def __init__(self, path):
        archive = tarfile.open(path)

        self.compounds = []
        self.embeddings = {}
        self.dict_atoms = []
        # Dataset reading

        for file in tqdm(archive.getmembers()):
            with archive.extractfile(file) as f:
                n_atoms = int(f.readline().decode("utf-8").strip())
                scalar = f.readline().decode("utf-8").split('\t')

                formula = {}

                for i in range(n_atoms):
                    atom = f.readline().decode("utf-8").split('\t')[0]

                    if atom not in self.dict_atoms:
                        self.dict_atoms.append(atom)

                    if atom in formula:
                        formula[atom] += 1
                    else:
                        formula[atom] = 1

                gibbs_energy = round(float(scalar[14]) * 627.503, 2)  # Перевод в kcal/mol
                smiles = f.readlines()[-2].decode("utf-8").strip().split('\t')

                self.compounds.append(
                    Compound(
                        formula=formula,
                        gibbs_energy=gibbs_energy,
                        smiles=smiles,
                        dataset_presence=set()))

        for atom in self.dict_atoms:
            for compound in self.compounds:
                if atom not in compound.formula.keys():
                    compound.formula[atom] = 0

    def boolean_map(
            self,
            dataset_path,
            dataset_name='pubchem'):
        """Checks each compound in the dataset for presence in another dataset
        """
        set_dataset = set()

        datasets = [dataset for dataset in os.listdir(dataset_path)]

        for dataset in tqdm(datasets):
            with gzip.open(os.path.join(dataset_path, dataset), 'rt') as f:
                for line_inchi in f:
                    one_inchi = line_inchi.strip('\n')
                    set_dataset.add(one_inchi)

        for compound in self.compounds:
            inchi = Chem.MolToInchi(compound.mol_structure)
            compound.dataset_presence.add(f'{dataset_name}_{inchi in set_dataset}')

        del set_dataset


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--ds-path', dest='dataset_path', type=str,
                        help='dataset path of compounds for generation reactions')
    parser.add_argument('--db-path', dest='database_path', type=str,
                        help='database path (for example, pubchem or zinc)')
    parser.add_argument('--db-name', dest='database_name', type=str,
                        help='database name (for example, "pubchem" or "zinc")')
    parser.add_argument('--n-jobs', dest='n_jobs', type=int,
                        help='Cores number for finder')
    parser.add_argument('--output-dir', dest='path_to_sep_dir_cmpds', type=str,
                        help='dir path of separated compounds')
    parser.add_argument('--output', dest='path_to_output_cmpds', type=str,
                        help='binary file path of treat compounds')
    args = parser.parse_args()

    kilomolecules = KiloMolecules(args.dataset_path)
    kilomolecules.boolean_map(args.database_path, args.database_name)

    split_cmpds = [list(sep_cmpds) for sep_cmpds in np.array_split(kilomolecules.compounds, args.n_jobs)]

    for i, sep_cmpds in enumerate(split_cmpds):
        with open(f"{os.path.join(args.path_to_sep_dir_cmpds, f'sep_cmpds_{i}.pkl')}", 'wb') as f:
            pkl.dump(sep_cmpds, f)

    with open(args.path_to_output_cmpds, 'wb') as f:
        pkl.dump(kilomolecules, f)
