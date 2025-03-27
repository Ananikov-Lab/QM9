from openbabel import pybel
from pybel import *


def getreversed(smiles):
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("smi")
    mol = pybel.readstring("smi", smiles)
    dimer = pybel.readstring("smi", smiles + smiles)

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


list_test_for_reverse = ['C=CCCCC1CC1C#C', 'C2CCc1ccccc1CC2', 'CCNC1(Br)CCC1ClCS']
list_some_test = ['c1ccccc1C#CCN', 'C#CCCCO', 'n1cccc2ccccc12', 'C=CCCCC#C',
                  'C\C=C/CC', 'C\\C=C\\CCCC#CC=C']
list_test_reverse = []
list_some_reversed = ['NCC#Cc1ccccc1', 'OCCCC#C', 'c12ccccc2cccn1', 'C#CCCCC=C',
                      'CC/C=C\C', 'C=CC#CCCC/C=C/C']
list_test_reverse_some = []

for cmpd in list_test_for_reverse:
    reverse_cmpd = getreversed(cmpd)
    list_test_reverse.append(reverse_cmpd)

for cmpd in list_some_test:
    reverse_cmpd = getreversed(cmpd)
    list_test_reverse_some.append(reverse_cmpd)


def test_cnt_presmiles_cmpd():
    for i in range(len(list_test_reverse)):
        mol = pybel.readstring('smi', list_test_for_reverse[i])
        smiles = mol.write('smi')
        assert len(list_test_reverse[i]) == len(smiles)


def test_cnt_presmiles_cmpd_c():
    uniq = {'C', 'N', 'S', '=', '#', 'n', '1', '2', 'c', 'B', 'r', 'l', ']', '[',
            '(', 'i', ')', '3'}
    for i in range(len(list_test_reverse)):
        for smiles in list_test_reverse[i]:
            assert (set(smiles) & uniq) == set(smiles)


def test_some_examples():
    for i in range(len(list_some_reversed)):
        mol = pybel.readstring('smi', list_some_reversed[i])
        inchi_noreverse = mol.write('inchi')
        mol_ = pybel.readstring('smi', list_test_reverse_some[i])
        inchi_test = mol_.write('inchi')
        print(inchi_test, inchi_noreverse)
        assert inchi_test == inchi_noreverse
