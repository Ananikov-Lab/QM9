from urllib.request import urlopen
from urllib.parse import quote
import re
from tqdm import tqdm
from rdkit import Chem
from argparse import ArgumentParser


def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'


def lab_alkynes(path_input, path_output):
    with open(path_input, 'rt', encoding='windows-1252') as f:
        acet = Chem.MolFromSmarts('C#C')
        lines = f.readlines()
        for i in tqdm(range(len(lines))):
            line = lines[i].split('\n')
            if len(re.findall('<Mol_ID>', line[0])) > 0:
                line_index = lines[i + 1].split('\n')
                index = line_index[0]
                line_cas = lines[i + 12].split('\n')
                if len(re.findall('<CAS', line_cas[0])) > 0:
                    line_cas_ = lines[i + 13].split('\n')
                    cas = line_cas_[0]
                    smiles = CIRconvert(cas)
                    if smiles == 'Did not work':
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    try:
                        if len(mol.GetSubstructMatches(acet)) > 0:
                            cmpd = index + '\t' + smiles + '\n'
                            with open(path_output, 'at') as f_db:
                                f_db.write(cmpd)
                    except Exception as err:
                        print(err)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input-db', dest='input', type=str, help='input of lab database')
    parser.add_argument('--output-db', dest='output', type=str, help='output of lab database')
    args = parser.parse_args()

    lab_alkynes(args.input, args.output)
