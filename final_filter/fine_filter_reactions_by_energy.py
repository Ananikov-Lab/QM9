from argparse import ArgumentParser
import pandas as pd


def filter_reactions(input_csv, n, output_csv):
    need_reactions = []
    df = pd.read_csv(input_csv, sep=';')
    df.sort_values(by=['gibbs energy'], inplace=True, ascending=True)
    df_need = df.iloc[:n]
    df_need.to_csv(output_csv)

    return need_reactions


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input-csv', dest='csv', type=str,
                        help='csv file with calculated reactions')
    parser.add_argument('--output', dest='output_need_reactions', type=str,
                        help='output path csv of need reactions for experiments')
    parser.add_argument('--number', dest='number', type=int,
                        help='number of need reactions')
    args = parser.parse_args()

    filter_reactions(args.csv, args.number, args.output_need_reactions)
