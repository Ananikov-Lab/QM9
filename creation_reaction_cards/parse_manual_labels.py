from argparse import ArgumentParser
import cv2
import os
import zipfile
import pandas as pd
import numpy as np
import pickle as pkl


def extract_zip(path_to_zip, output_zip):
    archive = zipfile.ZipFile(path_to_zip, 'r')
    archive.extractall(path=output_zip)


def marks_parser(output_numbers, closing):
    dict_outputs = {
        'id reaction': output_numbers,
        'Run': [],
        'Maybe run': [],
        'Not run': [],
        'Metal': [],
        'Acid': [],
        'Base': [],
        'No': [],
        'Soft': [],
        'Middle': [],
        'Hard': [],
        'In literature: yes': [],
        'In literature: no': [],
        'Synthetic utility': [],
        'The best reaction': []
    }
    for i in range(len(output_numbers)):
        if i == 0:
            m_x = 0
            m_y = 0
        if i == 1:
            m_x = 278
            m_y = 0
        if i == 2:
            m_x = 0
            m_y = 191
        if i == 3:
            m_x = 278
            m_y = 191
        set_react_run_ = set()
        set_react_mb_ = set()
        set_react_no_ = set()

        set_cat_met_ = set()
        set_cat_acid_ = set()
        set_cat_base_ = set()
        set_cat_nothing_ = set()

        set_condit_soft_ = set()
        set_condit_middle_ = set()
        set_condit_hard_ = set()

        set_yes_ = set()
        set_no_ = set()

        set_synthetic_ = set()

        set_best_rxn_ = set()
        if closing[98 + m_y, 43 + m_x] == 0:
            set_react_run_.add(1)
            set_react_mb_.add(0)
            set_react_no_.add(0)
        elif closing[113 + m_y, 43 + m_x] == 0:
            set_react_run_.add(0)
            set_react_mb_.add(1)
            set_react_no_.add(0)
        elif closing[129 + m_y, 43 + m_x] == 0:
            set_react_run_.add(0)
            set_react_mb_.add(0)
            set_react_no_.add(1)
        else:
            set_react_run_.add(0)
            set_react_mb_.add(0)
            set_react_no_.add(0)

        if closing[98 + m_y, 116 + m_x] == 0:
            set_cat_met_.add(1)
            set_cat_acid_.add(0)
            set_cat_base_.add(0)
            set_cat_nothing_.add(0)
        elif closing[113 + m_y, 116 + m_x] == 0:
            set_cat_met_.add(0)
            set_cat_acid_.add(1)
            set_cat_base_.add(0)
            set_cat_nothing_.add(0)
        elif closing[129 + m_y, 116 + m_x] == 0:
            set_cat_met_.add(0)
            set_cat_acid_.add(0)
            set_cat_base_.add(1)
            set_cat_nothing_.add(0)
        elif closing[144 + m_y, 116 + m_x] == 0:
            set_cat_met_.add(0)
            set_cat_acid_.add(0)
            set_cat_base_.add(0)
            set_cat_nothing_.add(1)
        else:
            set_cat_met_.add(0)
            set_cat_acid_.add(0)
            set_cat_base_.add(0)
            set_cat_nothing_.add(0)

        if closing[98 + m_y, 198 + m_x] == 0:
            set_condit_soft_.add(1)
            set_condit_middle_.add(0)
            set_condit_hard_.add(0)
        elif closing[113 + m_y, 198 + m_x] == 0:
            set_condit_soft_.add(0)
            set_condit_middle_.add(1)
            set_condit_hard_.add(0)
        elif closing[129 + m_y, 198 + m_x] == 0:
            set_condit_soft_.add(0)
            set_condit_middle_.add(0)
            set_condit_hard_.add(1)
        else:
            set_condit_soft_.add(0)
            set_condit_middle_.add(0)
            set_condit_hard_.add(0)

        if closing[163 + m_y, 182 + m_x] == 0:
            set_yes_.add(1)
            set_no_.add(0)
        elif closing[163 + m_y, 206 + m_x] == 0:
            set_yes_.add(0)
            set_no_.add(1)
        else:
            set_yes_.add(0)
            set_no_.add(0)

        if closing[182 + m_y, 126 + m_x] == 0:
            set_synthetic_.add(1)
        elif closing[182 + m_y, 149 + m_x] == 0:
            set_synthetic_.add(2)
        elif closing[182 + m_y, 174 + m_x] == 0:
            set_synthetic_.add(3)
        elif closing[182 + m_y, 199 + m_x] == 0:
            set_synthetic_.add(4)
        elif closing[182 + m_y, 223 + m_x] == 0:
            set_synthetic_.add(5)
        else:
            set_synthetic_.add(0)

        if closing[197 + m_y, 190 + m_x] == 0:
            set_best_rxn_.add(1)
        else:
            set_best_rxn_.add(0)

        dict_outputs['Run'].append(int(list(set_react_run_)[0]))
        dict_outputs['Maybe run'].append(int(list(set_react_mb_)[0]))
        dict_outputs['Not run'].append(int(list(set_react_no_)[0]))

        dict_outputs['Metal'].append(int(list(set_cat_met_)[0]))
        dict_outputs['Acid'].append(int(list(set_cat_acid_)[0]))
        dict_outputs['Base'].append(int(list(set_cat_base_)[0]))
        dict_outputs['No'].append(int(list(set_cat_nothing_)[0]))

        dict_outputs['Soft'].append(int(list(set_condit_soft_)[0]))
        dict_outputs['Middle'].append(int(list(set_condit_middle_)[0]))
        dict_outputs['Hard'].append(int(list(set_condit_hard_)[0]))

        dict_outputs['In literature: yes'].append(int(list(set_yes_)[0]))
        dict_outputs['In literature: no'].append(int(list(set_no_)[0]))

        dict_outputs['Synthetic utility'].append(int(list(set_synthetic_)[0]))

        dict_outputs['The best reaction'].append(int(list(set_best_rxn_)[0]))
    return dict_outputs


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--archive', dest='archive_labels', type=str,
                        help='zip of labels')
    parser.add_argument('--output', dest='output_labels', type=str,
                        help='output extract zip')
    parser.add_argument('--numbers', dest='output_numbers', type=str,
                        help='path of output numbers reaction id in zip')
    parser.add_argument('--csv', dest='csv_output', type=str,
                        help='path of csv output')
    args = parser.parse_args()

    extract_zip(args.archive_labels, args.output_labels)
    
    filenames = [filename for filename in os.listdir(args.output_labels)]
    
    for filename in filenames:
        i = filename.split('.')[-2][-2:]
        os.rename(os.path.join(args.output_labels, filename), os.path.join(args.output_labels, f'reactions-{i}.jpg'))
        
    i = 4
    with open(args.output_numbers, 'rb') as f:
        output_numbers = pkl.load(f)
    
    kernel = np.ones((7, 7), np.uint8)
    img = cv2.imread(f"{os.path.join(args.output_labels, 'reactions-01.jpg')}", cv2.COLOR_BGR2HSV)
    closing = cv2.morphologyEx(img[:, :, 0], cv2.MORPH_CLOSE, kernel)
    
    for x in range(img.shape[1]):
        for y in range(img.shape[0]):
            color = closing[y, x]
            if color < 240:
                closing[y, x] = 0
                    
    df = pd.DataFrame(marks_parser(output_numbers[:4], closing))

    for ind in range(2, len(os.listdir(args.output_labels)) + 1):
        print(ind)
        if ind < 10:
            ind = '0' + str(ind)
        img = cv2.imread(f"{os.path.join(args.output_labels, f'reactions-{ind}.jpg')}", cv2.COLOR_BGR2HSV)
        closing = cv2.morphologyEx(img[:, :, 0], cv2.MORPH_CLOSE, kernel)
        for x in range(img.shape[1]):
            for y in range(img.shape[0]):
                color = closing[y, x]
                if color < 240:
                    closing[y, x] = 0

        output_numbers_ = output_numbers[i:i + 4]
        df1 = pd.DataFrame(marks_parser(output_numbers_, closing))
        df = pd.concat([df, df1])
        i += 4
        
    df.reset_index(inplace=True, drop=True)
    df.to_csv(args.csv_output)
