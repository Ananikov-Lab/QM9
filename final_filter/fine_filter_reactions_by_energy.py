from argparse import ArgumentParser
import pickle as pkl
import pandas as pd
import numpy as np

def filter_reactions(input_csv):
    need_reactions = []
    df = pd.read_csv(input_csv)
    df.sort_values(by=['gibbs energy'], inplace=True, ascending=False).iloc[:10]
    
    for row in df.rows:
        need_reactions.append((row['reag'], row['alkyne'], row['prod']))

    return need_reactions



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input-csv', dest='csv', type=str, 
                       help='csv file with calculated reactions')
    parser.add_argument('--output', dest='output_need_reactions', type=str, 
                       help='output path of need reactions for experiments')
    args = parser.parse_args()
    
    reactions_for_expetiments = filter_reactions(args.input_csv)
   
    with open(args.output, 'wb') as f:
        pkl.dump(reactions_for_experiment, f)
        