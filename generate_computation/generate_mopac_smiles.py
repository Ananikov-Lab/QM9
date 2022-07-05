import os
from argparse import ArgumentParser

import numpy as np
from tqdm import tqdm

from openbabel import pybel
from pybel import *


def getreversed(smiles):
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("smi")
    mol = pybel.readstring("smi", smiles)
    dimer = pybel.readstring("smi", smiles + smiles)

    # Identify the 'last' atom (not necessarily the last
    # atom in the SMILES string)
    # One way to do this is to find what atom has a new
    # connection in the dimer
    last = 0
    for i, (a, b) in enumerate(zip(mol, dimer)):
        if a.OBAtom.GetExplicitDegree() != b.OBAtom.GetExplicitDegree():
            last = i + 1
            break
    assert last != 0

    # Generate the reverse
    conv.SetOptions('f"%d"l"%d"' % (last, 1), conv.OUTOPTIONS)
    reverse = conv.WriteString(mol.OBMol).rstrip()

    # Sanity check!! Dimer of 'reverse' should be identical
    newdimer = pybel.readstring("smi", reverse + reverse)
    assert dimer.write("can") == newdimer.write("can")

    return reverse


def parse_file_rxn(path):
    list_rxns = []
    with open(path, 'rt') as f:
        lines = f.readlines()
        for line in lines:
            line = line[:-1]
            line = line.replace(']1', ']5')
            line = line.replace(']2', ']6')
            list_rxns.append(line)

    return list_rxns


def parse_list_prods(list_rxn):
    list_prods_smiles = []
    for rxn in list_rxn:
        for i in range(len(rxn) - 1):
            if rxn[i] == '>':
                list_prods_smiles.append(rxn[i + 2:])
                break

    return list_prods_smiles


def generate_rxn_smiles(path_to_input, path_to_output, list_prods_smiles):
    set_smiles = set()
    with open(path_to_input, 'rt') as f:
        lines = f.readlines()
        for line in tqdm(lines):
            line_ = line.split('\t')
            index = line_[0]
            smiles = line_[1][:-1]
            list_smiles = list(smiles)
            bond_indexs = []
            for i in range(len(list_smiles) - 2):
                if list_smiles[i] == 'C' \
                        and list_smiles[i + 1] == '#' \
                        and list_smiles[i + 2] == 'C':
                    bond_indexs += [i]
            for indexs in bond_indexs:
                r1 = ''
                r2 = ''
                bond_index = indexs
                ind = bond_indexs.index(indexs)
                try:
                    r1 = f'({smiles[:bond_index]})'
                except:
                    r1 = 'H'
                    print(1)

                try:
                    r2 = f'({smiles[bond_index + 3:]})'
                except:
                    r2 = 'H'

                if r1 == '()' or r1 == '':
                    r1 = '([H])'
                if r1 == '([H])':
                    r1_mod = '([H])'
                else:
                    r1_mod = getreversed(r1[1:-1])
                if r2 == '()' or r2 == '':
                    r2 = '([H])'

                for cmpd in list_prods_smiles:
                    m = list_prods_smiles.index(cmpd)
                    cmpd = cmpd.replace('([Cl:1])', f'({r1_mod})')
                    cmpd = cmpd.replace('[Cl:1]', f'({r1_mod})')
                    cmpd = cmpd.replace('([Cl:4])', f'{r2}')
                    cmpd = cmpd.replace('[Cl:4]', f'{r2}')
                    mol = pybel.readstring('smi', cmpd)
                    inchi = mol.write('inchi')
                    mol = pybel.readstring('inchi', inchi)
                    cmpd = mol.write('smi')[:-2]
                    if cmpd not in set_smiles:
                        set_smiles.add(cmpd)
                        with open(path_to_output, 'at') as fil:
                            cmpd = f'{m}_' + index + f'_{ind}' + '_1\t' + cmpd + '\n'
                            fil.write(cmpd)

                for cmpd in list_prods_smiles:
                    m = list_prods_smiles.index(cmpd)
                    cmpd = cmpd.replace('([Cl:1])', f'{r2}')
                    cmpd = cmpd.replace('[Cl:1]', f'{r2}')
                    cmpd = cmpd.replace('([Cl:4])', f'({r1_mod})')
                    cmpd = cmpd.replace('[Cl:4]', f'({r1_mod})')
                    mol = pybel.readstring('smi', cmpd)
                    inchi = mol.write('inchi')
                    mol = pybel.readstring('inchi', inchi)
                    cmpd = mol.write('smi')[:-2]
                    if cmpd not in set_smiles:
                        set_smiles.add(cmpd)
                        with open(path_to_output, 'at') as fil:
                            cmpd = f'{m}_' + index + f'_{ind}' + '_2\t' + cmpd + '\n'
                            fil.write(cmpd)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('path_map_rxns', type=str,
                        help="First path of file with list mapper rxns")
    parser.add_argument('path_to_input', type=str,
                        help="Second - path of file with alkyne`s smiles")
    parser.add_argument('path_to_output', type=str,
                        help="Third - path of output file")
    args = parser.parse_args()

    list_rxns = parse_file_rxn(args.path_map_rxns)
    list_prods_smiles = parse_list_prods(list_rxns)
    generate_rxn_smiles(args.path_to_input, args.path_to_output, list_prods_smiles)
    