from find_reactions import *
from treatment_dataset import *

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--n-jobs', dest='n_jobs', type=int,
                        help='Cores number for finder')
    parser.add_argument('--output', dest='reactions', type=str,
                        help='Path to final reactions list')
    args = parser.parse_args()
    
    reactions = []
    
    for i in range(args.n_jobs):
        with open(f'reactions_list_{i}.pkl', 'rb') as f:
            reaction = pkl.load(f)
        reactions += reaction
    
    with open(args.reactions, 'wb') as f:
        pkl.dump(reactions, f)
        