from openbabel import pybel
from pybel import *
import os
from argparse import ArgumentParser
from tqdm import tqdm


def gjf_generate(path_to_input_dir, path_to_chk, ram, basis,
                 n_proc, path_to_output_dir):
    filenames = [filename for filename in os.listdir(path_to_input_dir) if filename.endswith('.out')]
    for filename in tqdm(filenames):
        with open(f'{path_to_input_dir}/{filename}') as myfile:
            lines = myfile.readlines()
            index = filename[:-4]
            with open(f'{path_to_output_dir}/multiprocess_short_{index}.gjf', 'wt') as f:
                f.write(
                    f'%chk={path_to_chk}/multiprocess_short_{index}.chk\n%NProcShared={n_proc}\n%Mem={ram}Mb\n# {basis} \n\nTitle Card Required\n\n0 1\n')
                for i in range((len(lines) - 1), 0, -1):
                    if lines[i].startswith('           Empirical Formula:'):
                        i -= 3
                        for m in range(i, 0, -1):
                            f.write(lines[m][7:])
                            if lines[m].startswith('   1    '):
                                break
                        break


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('path_to_input_dir', type=str,
                        help="First - path of dir with mopac`s outputs")
    parser.add_argument('path_to_chk', type=str,
                        help="Second - path of chk")
    parser.add_argument('ram', type=str,
                        help="Third - ram")
    parser.add_argument('basis', type=str,
                        help="Fourth - basis of calculation")
    parser.add_argument('n_proc', type=str,
                        help="Fifth - No of procs")
    parser.add_argument('path_to_output_dir', type=str,
                        help="Sixth - path of output dir")
    args = parser.parse_args()

    gjf_generate(args.path_to_input_dir, args.path_to_chk,
                 args.ram, args,basis, args.n_proc, args.path_to_output_dir)
    