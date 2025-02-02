from argparse import ArgumentParser
from tqdm import tqdm
import os

def energy(path_to_mols, smiles_prods_dict, output_file):
    filenames = [filename for filename in os.listdir(path_to_mols) if filename.endswith('.log')]
    for filename in tqdm(filenames):
        index = filename.split('_', 2)[-1]
        with open(os.path.join(path_to_mols, filename)) as myfile:
            lines = myfile.readlines()
            for line in lines:
                if line.startswith(' Sum of electronic and thermal Free Energies= '):
                    line = line.strip().split(' ')
                    mol_energy = float(line[-1])
                    smiles = smiles_prods_dict[index]
                    with open(output_file, 'at') as f:
                        f.write(f'{index}\t{smiles}\t{mol_energy}\n')
            
            
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('path_to_mols', type=str,
                        help="First - path of dir file with molecules data")
    parser.add_argument('prods_smiles', type=str,
                        help="Second - path of .txt file with prods smiles")
    parser.add_argument('output_file', type=str,
                        help="Third - path of output .txt file with mols energy")
    args = parser.parse_args()
    
    smiles_dict = {}
    
    with open(args.prods_smiles, 'rt') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            index = line[0]
            smiles = line[1]
            smiles_dict[index] = smiles
    
    energy(args.path_to_mols, smiles_dict, args.output_file)
