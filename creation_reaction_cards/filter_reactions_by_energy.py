from argparse import ArgumentParser
import pickle as pkl

def filter_reactions(reactions, model):
    need_reactions = []
    label = model.labels_
    reaction_list, label = zip(*sorted(zip(reactions.mpairs_of_reacts, label)))
    flag = 0
    for i in range(len(reaction_list)):
        if reaction_list[i].hydro_shift == 'Without hydro-shift':
            need_reactions.append(reaction_list[i])
            continue
        if m == 18:
            flag = 0
            continue
        if reaction_list[i].reac_rating == 'High probably':
            need_reactions.append(reaction_list[i])
            flag += 1
    return need_reactions

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--reactions', dest='input_reactions', type=str, 
                       help='binary file of generated reactions')
    parser.add_argument('--model', dest='model', type=str, 
                       help='model`s file')
    parser.add_argument('--output', dest='output_need_reactions', type=str, 
                       help='output path of need reactions for expert verification')
    args = parser.parse_args()
    
    with open(args.input_reactions, 'rb') as f:
        reactions = pkl.load(f)
    
    with open(args.model, 'rb') as f:
        sc_model = pkl.load(f)
    
    need_reactions_list = filter_reactions(reactions, sc_model)
    
    with open(args.output_need_reactions, 'wb') as f:
        pkl.dump(need_reactions_list, f)
        