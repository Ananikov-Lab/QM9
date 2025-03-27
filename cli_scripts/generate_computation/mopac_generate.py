from openbabel import pybel
from pybel import *
from tqdm import tqdm
from argparse import ArgumentParser


def generate_mopac(path_to_input, path_to_output_dir):
    with open(path_to_input, 'rt') as f:
        lines = f.readlines()
        for line in tqdm(lines):
            line_ = line.split('\t')
            index = line_[0]
            smiles = line_[1][:-1]
            mol = pybel.readstring('smi', smiles)
            mol.make3D()
            coordinates = mol.write('gjf').split('\n')
            with open(f'{path_to_output_dir}/{index}.mop', 'wt') as fil:
                fil.write(
                    f'PM7 OPT BONDS AUX GRAPHF GNORM=0.01 + \nCHARGE=0\n {index}\n Optimize geometry and save extra info\n')
                for string in coordinates[6:-2]:
                    string_ = string + '\n'
                    fil.write(string_)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('path_to_input', type=str,
                        help="First - path of file with prods` smiles")
    parser.add_argument('path_to_output_dir', type=str,
                        help="Second - path of output dir")
    args = parser.parse_args()

    generate_mopac(args.path_to_input, args.path_to_output_dir)
