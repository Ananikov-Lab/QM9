from processing_dataset import *


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

    def plot(self):
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

    def generate_reactions(self, dataset, list_smarts, output_path_reacts, carb_cnt, alk_index, cmpds_list):
        substrate = dataset.compounds[alk_index]
        for i, compound in enumerate(dataset.compounds):
            self.compounds_id_dict[compound.inchi] = i
            
        set_unique_react = set()
        for compound_reag in tqdm(cmpds_list):
            if compound_reag.formula['C'] <= carb_cnt:
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
                                except KekulizeException:
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
            self.other_prp['Formation/abolition of functional group(s)'] = {}

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
                    self.other_prp['Formation/abolition of heterocycle(s)'][
                        name_hetero] = 'Heterocycle(s) is(are) not formed/abolished in react'
                elif difference_hetero < 0:
                    self.other_prp['Formation/abolition of heterocycle(s)'][
                        name_hetero] = 'Heterocycle(s) is(are) formed in react'
                elif difference_hetero > 0:
                    self.other_prp['Formation/abolition of heterocycle(s)'][
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


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--sep-file', dest='sep_file', type=str,
                        help='file path of compounds lists for generation reactions')
    parser.add_argument('--ds-path', dest='dataset_path', type=str,
                        help='dataset path of compounds for generation reactions')
    parser.add_argument('--carb', dest='carb_cnt', type=int,
                        default=7, help='max count of carbons in cmpd')
    parser.add_argument('--index', dest='alk_index', type=int,
                        default=3, help='index of desired alkyne in dataset')
    parser.add_argument('--ls-path', dest='list_smarts_path', type=str,
                        help='file path with reaction templates')
    parser.add_argument('--ir-path', dest='path_id_reacts', type=str,
                        help='file path with generated reactions using id compounds')
    parser.add_argument('--output', dest='path_to_output_list_reacts', type=str,
                        help='binary file path of generated reactions')
    args = parser.parse_args()

    with open(args.list_smarts_path, 'rb') as f:
        list_smarts = pkl.load(f)

    with open(args.dataset_path, 'rb') as f:
        kilomolecules = pkl.load(f)
     
    with open(args.sep_file, 'rb') as f:
        sep_cmpds = pkl.load(f)
    
    reactions = ReactionSearcher(kilomolecules)
    reactions.generate_reactions(kilomolecules, list_smarts, args.path_id_reacts, args.carb_cnt, args.alk_index,
                                 sep_cmpds)
    reactions.list_maker_reacts(args.path_id_reacts, kilomolecules)
    reactions.classifier()
    with open(args.path_to_output_list_reacts, 'wb') as f:
        pkl.dump(reactions, f)
