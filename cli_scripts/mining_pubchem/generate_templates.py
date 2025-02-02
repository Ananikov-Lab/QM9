from argparse import ArgumentParser
import pickle as pkl
from rdkit import Chem
from rdkit.Chem import AllChem


def cmpd_generate(smiles, cnt, type_cmpd='non-cycle'):
    """There is an enumeration of all possible options for the multiplicity of bonds within the compounds"""
    list_smiles = []
    if cnt == 0:
        if type_cmpd == 'non-cycle':
            cmpd = smiles + '*'
            list_smiles += [cmpd]
            return list_smiles
        elif type_cmpd == 'cycle':
            for end_bond in ['#', '=', '']:
                cmpd = '*1' + smiles[1:] + '*' + end_bond + '1'
                list_smiles += [cmpd]
            return list_smiles
        else:
            raise ValueError
    else:
        for bond in ['#', '=', '']:
            cmpd = smiles + '*' + bond
            list_smiles += cmpd_generate(cmpd, cnt - 1, type_cmpd)
            cmpd = smiles
        return list_smiles


def screening_cmpd(list_smiles):
    """Filtering the resulting layouts from the impossible multiplicity of atoms"""
    scn_smiles = []
    for smiles in list_smiles:
        if smiles.rfind('#*=') != -1 or smiles.rfind('#*#') != -1 \
                or smiles.rfind('=*#') != -1 or (smiles.rfind('*=1') != -1 and smiles.find('*1#') != -1) \
                or (smiles.rfind('*#1') != -1 and smiles.find('*1=') != -1) \
                or (smiles.rfind('*#1') != -1 and smiles.find('*1#') != -1):
            continue
        else:
            scn_smiles += [smiles]
    return scn_smiles


def atom_valent(list_smiles, type_cmpd='non-cycle'):
    """Creation of an additional list that will store the valencies of each of the atoms
    (needed to account for thehydride shift) """
    valent_bond = {
        '#': 3,
        '=': 2,
        '*': 1,
        '1': 1,
        '0': 0,
    }
    list_valent = []
    for smiles in list_smiles:
        val_smiles = []
        if type_cmpd == 'non-cycle':
            smiles_ = '0' + smiles + '0'
        elif type_cmpd == 'cycle':
            smiles_ = '1' + smiles.replace('1', '') + '1'
        else:
            raise ValueError
        for i in range(len(smiles_) - 1):
            if smiles_[i] == '*':
                val_smiles.append(valent_bond[smiles_[i - 1]] + valent_bond[smiles_[i + 1]])
        list_valent += [val_smiles]
    return list_valent


def atom_valent_rxn(list_valent):
    """Creating a list of valences inside a reaction"""
    atoms_valent = []
    for valent_ in list_valent:
        atoms_valent += [valent_ + [3, 3]]  # Acetylene
    return atoms_valent


def atom_numerate(list_smiles):
    """Creating smiles with numbered atoms"""
    list_numsmiles = []
    for smiles in list_smiles:
        numsmiles = ''
        i = 1
        for char in smiles:
            if char == '*':
                numsmiles += f'[{char}:{i}]'
                i += 1
            else:
                numsmiles += char
        list_numsmiles += [numsmiles]
    return list_numsmiles


def rxn_generate(list_numsmiles, list_numsmiles_c, valent_reag_w_acet, valent_prod):
    """Generation of mock reactions by "gluing" reactants and products.
    As well as taking into account the hydride shift"""
    list_rxn_smarts = []
    for smarts_reag, val_reag in zip(list_numsmiles, valent_reag_w_acet):
        i = int(smarts_reag[-2])
        smarts_frag = Chem.MolFromSmarts(smarts_reag)
        for smarts_prod, val_prod in zip(list_numsmiles_c, valent_prod):
            rxn_smarts = f'{smarts_reag}.[#6:{i + 1}]#[#6:{i + 2}]>>{smarts_prod}'
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            if val_reag == val_prod:
                hydro_shift = 'Without hydro-shift'
            else:
                hydro_shift = 'With hydro-shift'
            list_rxn_smarts += [[rxn, hydro_shift, smarts_frag]]
    return list_rxn_smarts


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--output', dest='output_path', type=str, help='path of output file')
    args = parser.parse_args()

    smiles = ''
    list_smarts = []
    for i in range(6):
        all_smiles = cmpd_generate(smiles, i, type_cmpd='non-cycle')
        all_smiles_c = cmpd_generate(smiles, i + 2, type_cmpd='cycle')

        list_normsmiles = screening_cmpd(all_smiles)
        list_normsmiles_c = screening_cmpd(all_smiles_c)

        valenting_list = atom_valent(list_normsmiles, type_cmpd='non-cycle')
        valenting_list_c = atom_valent(list_normsmiles_c, type_cmpd='cycle')
        valenting_list_w_acet = atom_valent_rxn(valenting_list)

        list_numsmiles_c = atom_numerate(list_normsmiles_c)
        list_numsmiles = atom_numerate(list_normsmiles)

        list_smarts.append(rxn_generate(list_numsmiles, list_numsmiles_c, valenting_list_w_acet, valenting_list_c))

    with open(args.output_path, 'wb') as f:
        pkl.dump(list_smarts, f)
