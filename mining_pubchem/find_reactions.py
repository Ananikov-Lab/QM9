import tarfile
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import gzip
import os
import re
# import rdkit components
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from joblib import Parallel, delayed
import gc
from multiprocessing import Pool


# use IPythonConsole for pretty drawings
from rdkit.Chem.Draw import IPythonConsole

## The next line is commented out
### because GitHub does not render svg's embedded in notebooks
# IPythonConsole.ipython_useSVG=True
IPythonConsole.ipython_useSVG = False

# for flattening tuples and lists
from itertools import chain
from rdkit import RDLogger
import random
from argparse import ArgumentParser
import pickle as pkl


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


class Reaction:
    def __init__(self, reag, prod, dif, reac_energy, hydro_shift):
        self.reag = reag
        self.prod = prod
        self.dif = dif
        self.reac_energy = reac_energy
        self.cmpd_availability = set()
        self.hydro_shift = hydro_shift
        self.other_prp = {}

        if self.reac_energy > 10:
            self.reac_rating = 'Unlikely'
        elif 10 >= self.reac_energy > 0:
            self.reac_rating = 'May take place'
        elif 0 >= self.reac_energy > -10:
            self.reac_rating = 'Probably'
        elif self.reac_energy <= -10:
            self.reac_rating = 'High probably'

    def availability(self, dataset_name):
        if f'{dataset_name}_True' in self.reag[0].dataset_presence \
                and f'{dataset_name}_False' in self.prod[0].dataset_presence:
            self.cmpd_availability.add(f'{dataset_name}_only_reag')
        elif f'{dataset_name}_False' in self.reag[0].dataset_presence \
                and f'{dataset_name}_True' in self.prod[0].dataset_presence:
            self.cmpd_availability.add(f'{dataset_name}_only_prod')
        elif f'{dataset_name}_True' in self.reag[0].dataset_presence \
                and f'{dataset_name}_True' in self.prod[0].dataset_presence:
            self.cmpd_availability.add(f'{dataset_name}_reag+prod')
        else:
            self.cmpd_availability.add(f'{dataset_name}_both_miss')

    def plot(self, type_='First'):
        if type_ != 'First':
            raise NotImplementedError

        reagent = self.reag[0]
        product = self.prod[0]
        substrate = self.reag[1]

        print('React\'s Energy:', self.reac_energy, 'kcal/mol')
        print('React\'s Rating:', self.reac_rating)
        print('Availability:', self.cmpd_availability)
        print('Hydro-shift:', self.hydro_shift)

        rxn = rdChemReactions.ReactionFromSmarts(f'{reagent.smiles[0]}>{substrate.smiles[0]}>{product.smiles[0]}',
                                                 useSmiles=True)
        img = Chem.Draw.ReactionToImage(
            rxn, subImgSize=(300, 300)
        )
        plt.imshow(img)
        plt.axis('off')
        plt.show()


class ReactionSearcher:
    def __init__(self, dataset, dataset_name=None):
        self.compounds = []
        self.dataset_name = dataset_name
        self.pairs_of_reacts = []
        self.mpairs_of_reacts = []
        self.other_prp = {}
        self.compounds_id_dict = {}

        if self.dataset_name is not None:
            for compound in dataset.compounds:
                if f'{dataset_name}_True' in compound.dataset_presence:
                    self.compounds.append(compound)
        else:
            self.compounds = dataset.compounds

    def generate_reactions(self, dataset, list_smarts, output_path_reacts):
        substrate = dataset.compounds[3]  # Acetylene
        for i, compound in enumerate(dataset.compounds):
            self.compounds_id_dict[compound.inchi] = i

        set_unique_react = set()
        for compound_reag in tqdm(dataset.compounds):
            if compound_reag.formula['C'] <= 7:
                for rxn_number in list_smarts:
                    for rxn in rxn_number:
                        Chem.Kekulize(compound_reag.mol_structure)
                        if len(compound_reag.mol_structure.GetSubstructMatches(rxn[2])) > 0:
                            run_react = rxn[0].RunReactants(
                                (compound_reag.mol_structure, substrate.mol_structure))
                            for react in run_react:
                                try:
                                    react_smiles = Chem.MolToSmiles(react[0])
                                    react_norm_mol = Chem.MolFromSmiles(react_smiles)
                                    if react_norm_mol is not None:
                                        react_inchi = Chem.MolToInchi(react_norm_mol)
                                        if react_inchi in self.compounds_id_dict.keys():
                                            id_reag = self.compounds_id_dict[compound_reag.inchi]
                                            id_prod = self.compounds_id_dict[react_inchi]
                                            id_sub = self.compounds_id_dict[substrate.inchi]
                                            pair_of_reag_prod = f'{id_reag}\t{id_sub}\t{id_prod}\t{rxn[1]}\n'
                                            compound_prod = dataset.compounds[id_prod]
                                            H1, H2, H3 = compound_reag.formula['H'], \
                                                         compound_prod.formula['H'], substrate.formula['H']
                                            if (H2 - H1 == H3) and (pair_of_reag_prod not in set_unique_react):
                                                with open(output_path_reacts,
                                                          'at') as f:
                                                    f.write(pair_of_reag_prod)
                                                set_unique_react.add(pair_of_reag_prod)
                                except:
                                    pass

    def list_maker_reacts(self, path, dataset):
        list_all_reacts = []
        with open(path, 'rt') as f:
            for line in f:
                line_ = line.rstrip('\n')
                react_info = line_.split('\t')
                reag = dataset.compounds[int(react_info[0])]
                sub = dataset.compounds[int(react_info[1])]
                prod = dataset.compounds[int(react_info[2])]
                reac_energy = prod.gibbs_energy - reag.gibbs_energy - sub.gibbs_energy
                reaction = Reaction([reag, sub],
                                    [prod],
                                    1,
                                    reac_energy,
                                    react_info[3]  # Info of hydro-shift
                                    )
                reaction.availability('pubchem')
                reaction.availability('CAS')
                list_all_reacts.append(reaction)

        sorts_of_reacts = sorted(list_all_reacts,
                                 key=lambda react: react.reac_energy)

        reag = sorts_of_reacts[0].reag[0]
        prod = sorts_of_reacts[0].prod[0]
        hydro_shift = sorts_of_reacts[0].hydro_shift
        self.mpairs_of_reacts.append(sorts_of_reacts[0])
        for react in sorts_of_reacts:
            if reag == react.reag[0] and prod == react.prod[0] and hydro_shift == react.hydro_shift:
                continue
            else:
                self.mpairs_of_reacts.append(react)
                reag = react.reag[0]
                prod = react.prod[0]
                hydro_shift = react.hydro_shift
        del list_all_reacts, sorts_of_reacts

    def outputs_reacts(self):
        for react in self.mpairs_of_reacts:
            react.plot()

    def classifier(self):
        heretocycles = {
            'pyrrole': Chem.MolFromSmiles('C1=CNC=C1'),
            'furane': Chem.MolFromSmiles('C1=COC=C1'),
            'thiophene': Chem.MolFromSmiles('C1=CSC=C1')
        }

        all_func_group = {
            'NR3': Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]'),
            'CN': Chem.MolFromSmarts('[NX1]#[CX2]'),
            'COOH': Chem.MolFromSmarts('[CX3](=O)[OX2H1]'),
            'ROR': Chem.MolFromSmarts('[*][O][*]'),
            'OH': Chem.MolFromSmarts('[OX2H]'),
            'C(O)': Chem.MolFromSmarts('[CX3]=[OX1]')
        }

        for react in self.mpairs_of_reacts:
            reag = react.reag[0].mol_structure
            prod = react.prod[0].mol_structure
            self.other_prp['Formation/abolition of heterocycle(s)'] = {}
            self.other_prp['Formation/abolition of functional group'] = {}

            for i in range(10):
                if i == 0 and (Chem.GetSSSR(reag) - Chem.GetSSSR(prod)) == 0:
                    self.other_prp['Cycles in react'] = 'New cycle(s) is not formed'
                elif (Chem.GetSSSR(reag) - Chem.GetSSSR(prod)) == i:
                    self.other_prp['Cycles in react'] = f'{i} cycle(s) is formed'
                elif (Chem.GetSSSR(reag) - Chem.GetSSSR(prod)) == -i:
                    self.other_prp['Cycles in react'] = f'{i} cycle(s) is abolished'

            for name_hetero, mol_hetero in heretocycles.items():
                difference_hetero = len(reag.GetSubstructMatches(mol_hetero)) - len(
                    prod.GetSubstructMatches(mol_hetero))
                if difference_hetero == 0:
                    self.other_prp['Formation of heterocycle(s)'][
                        name_hetero] = 'Heterocycle(s) is(are) not formed/abolished in react'
                elif difference_hetero < 0:
                    self.other_prp['Formation of heterocycle(s)'][
                        name_hetero] = 'Heterocycle(s) is(are) formed in react'
                elif difference_hetero > 0:
                    self.other_prp['Formation of heterocycle(s)'][
                        name_hetero] = 'Heterocycle(s) is(are) abolished in react'

            for name_group, smarts_group in all_func_group.items():
                difference_group = len(reag.GetSubstructMatches(smarts_group)) - len(
                    prod.GetSubstructMatches(smarts_group))
                if difference_group == 0:
                    self.other_prp['Formation/abolition of functional group(s)'][
                        name_group] = 'Functional group(s) is(are) not formed/abolished in react'
                elif difference_group < 0:
                    self.other_prp['Formation/abolition of functional group(s)'][
                        name_group] = 'Functional group(s) is(are) formed in react'
                elif difference_group > 0:
                    self.other_prp['Formation/abolition of functional group(s)'][
                        name_group] = 'Functional group(s) is(are) abolished in react'


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

                gibbs_energy = round(float(scalar[14]) / 627.503, 2)  # Перевод в kcal/mol
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
            dataset_name='pubchem',
            mapping_mode='structure'):
        """Checks each compound in the dataset for presence in another dataset
        """
        if mapping_mode != 'structure':
            raise NotImplementedError

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
    parser.add_argument('dataset_path', type=str)
    parser.add_argument('database_path', type=str)
    parser.add_argument('database_name', type=str)
    parser.add_argument('list_smarts_path', type=str)
    parser.add_argument('path_id_reacts', type=str)
    parser.add_argument('path_to_output_list_reacts', type=str)
    args = parser.parse_args()

    kilomolecules = KiloMolecules(args.dataset_path)
    kilomolecules.boolean_map(args.database_path, args.database_name)

    with open(args.list_smarts_path, 'rb') as f:
        list_smarts = pkl.load(f)

    reactions = ReactionSearcher(kilomolecules)
    reactions.generate_reactions(kilomolecules, list_smarts, args.path_id_reacts)
    reactions.list_maker_reacts(args.path_id_reacts, kilomolecules)
    reactions.classifier()
    with open(args.path_to_output_list_reacts, 'wb') as f:
        pkl.dump(reactions, f)
