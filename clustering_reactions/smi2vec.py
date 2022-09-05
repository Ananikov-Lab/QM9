import torch
from torch import nn
from tqdm.notebook import tqdm
import numpy as np
from argparse import ArgumentParser
import sys
sys.path.append('../../smiles_transformer/smiles_transformer')

from pretrain_trfm import TrfmSeq2seq
from pretrain_rnn import RNNSeq2Seq
from build_vocab import WordVocab
from utils import split



def get_inputs(sm):
    seq_len = 220
    sm = sm.split()
    if len(sm)>218:
        print('SMILES is too long ({:d})'.format(len(sm)))
        sm = sm[:109]+sm[-109:]
    ids = [vocab.stoi.get(token, unk_index) for token in sm]
    ids = [sos_index] + ids + [eos_index]
    seg = [1]*len(ids)
    padding = [pad_index]*(seq_len - len(ids))
    ids.extend(padding), seg.extend(padding)
    return ids, seg

def get_array(smiles):
    x_id, x_seg = [], []
    for sm in smiles:
        a,b = get_inputs(sm)
        x_id.append(a)
        x_seg.append(b)
    return torch.tensor(x_id), torch.tensor(x_seg)

def process_smiles(sm):
    sm = ' '.join(list(sm))
    before = ['C l -', 'C l', 'O -', 'N +', 'n +', 'B r -', 'B r', 'N a +', 'N a', 'I -', 'S i']
    after = ['Cl-', 'Cl', 'O-', 'N+', 'n+', 'Br-', 'Br', 'Na+', 'Na', 'I-', 'Si']
    for b, a in zip(before, after):
        sm = sm.replace(b, a)
    
    return sm

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--p-vocab', dest='vocab', type=str, help='vocab file of transformer')
    parser.add_argument('--p-trfm', dest='trfm', type=str, help='path to transformer')
    parser.add_argument('--smi', dest='smiles', type=str, help='file with smiles')
    parser.add_argument('--output', dest='output_enb', type=str, help='npy output file')
    args = parser.parse_args()
    
    
    vocab = WordVocab.load_vocab(args.vocab)
    trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    trfm.load_state_dict(torch.load(args.trfm))
    trfm.eval()
    
    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    mask_index = 4
    smiles = []
    
    with open(args.smiles, 'rt') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            smi = line
            smiles.append(smi)
            
    xid, _ = get_array(smiles)
    X = trfm.encode(torch.t(xid))
    
    with open(args.output_enb, 'wb') as f:
        np.save(f, X)
