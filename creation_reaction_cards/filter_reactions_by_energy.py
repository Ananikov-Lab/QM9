import pickle as pkl
from mining_pubchem.find_reactions import *
import random

def filter_reactions(reactions, model, n):
    need_reactions = []
    label = model.labels_
    reaction_list, label = zip(*sorted(zip(reactions.mpairs_of_reacts, label), reverse=True, key=lambda a: a[1]))
    flag = 0
    label_ = label[0]
    for i in range(len(reaction_list)):
        if reaction_list[i].hydro_shift == 'Without hydro-shift':
            need_reactions.append(reaction_list[i])
            continue
        if flag == n:
            flag = 0
            label_ -= 1
            continue
        if reaction_list[i].reac_rating == 'High probably' and label_ == label[i]:
            need_reactions.append(reaction_list[i])
            flag += 1
    return need_reactions


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--reactions', dest='input_reactions', type=str, 
                       help='binary file of generated reactions')
    parser.add_argument('--model', dest='model', type=str, 
                       help='model`s file')
    parser.add_argument('--number', dest='number_clust_reacts', type=int, 
                       help='number of reacts in one cluster which need you')
    parser.add_argument('--output', dest='output_need_reactions', type=str, 
                       help='output path of need reactions for expert verification')
    parser.add_argument('--output-numbers', dest='output_shuffle_reactions', type=str, 
                       help='output path of index reactions for expert verification')
    args = parser.parse_args()
    
    with open(args.input_reactions, 'rb') as f:
        reactions = pkl.load(f)
    
    with open(args.model, 'rb') as f:
        sc_model = pkl.load(f)
    
    need_reactions_list = filter_reactions(reactions, sc_model, args.number_clust_reacts)
    
    dict_w_o_shuffle = {}
    dict_w_shuffle = {}
    
    output_numbers = []
    
    for rxn in need_reactions_list:
        dict_w_o_shuffle[str(need_reactions_list.index(rxn))] = rxn
    
    random.shuffle(need_reactions_list)
    
    for rxn in need_reactions_list:
        dict_w_shuffle[str(need_reactions_list.index(rxn))] = rxn
        
    for name_wo, val_wo in dict_w_o_shuffle.items():
        for name_w, val_w in dict_w_shuffle.items():
            if val_wo == val_w:
                output_numbers.append(name_w)
    
    with open(args.output_need_reactions, 'wb') as f:
        pkl.dump(need_reactions_list, f)
    
    with open(args.output_shuffle_reactions, 'wb') as f:
        pkl.dump(output_numbers, f)
        