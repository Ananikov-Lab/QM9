from find_reactions import *
from processing_dataset import *
import os


def download_smiles(smiles_list, smiles_path):
    with open(smiles_path, 'wt') as f:
        for cmpd in smiles_list:
                text = cmpd + '\n'
                f.write(text)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', dest='w_dir', type=str,
                        help='Path to work dir')
    parser.add_argument('--output-rs', dest='r_smiles', type=str,
                        help='Path to output smiles reagents')
    parser.add_argument('--output-ps', dest='p_smiles', type=str,
                        help='Path to output smiles products')
    parser.add_argument('--output', dest='reactions', type=str,
                        help='Path to final reactions list')
    args = parser.parse_args()
    
    reactions = []
    filenames = [filename for filename in os.listdir(args.w_dir) if filename.endswith('.pkl')]
    
    for i in range(len(filenames)):
        with open(f"{os.path.join(args.w_dir, f'reactions_list_{i}.pkl')}",
                  'rb') as f:
            reaction = pkl.load(f)
        reactions += reaction
        
    reags = set()
    prods = set()
    for rxn in reactions:
        reag = Chem.MolToSmiles(rxn.reag[0].mol_structure)
        prod = Chem.MolToSmiles(rxn.prod[0].mol_structure)
        reags.add(reag)
        prods.add(prod)
        
    reags = list(reags)
    prods = list(prods)
    
    
    
    for path, smiles in zip([args.r_smiles, args.p_smiles],
                            [reags, prods]):
        download_smiles(smiles, path)
            
    
    with open(args.reactions, 'wb') as f:
        pkl.dump(reactions, f)
        