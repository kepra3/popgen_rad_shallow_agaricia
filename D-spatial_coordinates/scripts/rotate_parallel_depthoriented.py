#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 31/5/22
@description: Rotate each plot parallel
"""

import numpy as np
import pandas as pd

PATH = '../results'

def import_coords_compile(PATH):
    WP05 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_kal_05m_20200214.csv'.format(PATH))
    WP10 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_kal_10m_20200214.csv'.format(PATH))
    WP20 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_kal_20m_20200214.csv'.format(PATH))
    CA05 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_cas_05m_20201212.csv'.format(PATH))
    CA10 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_cas_10m_20201212.csv'.format(PATH))
    CA20 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_cas_20m_20201212.csv'.format(PATH))
    SB05 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_sna_05m_20200303.csv'.format(PATH))
    SB10 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_sna_10m_20200303.csv'.format(PATH))
    SB20 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_sna_20m_20200303.csv'.format(PATH))
    SQ12 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_seb_10m_20201210.csv'.format(PATH))
    SQ20 = pd.read_csv('{}/rotated_scaled_annotations_X_cur_seb_20m_20201210.csv'.format(PATH))

    all_plus90 = pd.concat([SB20, SQ20])
    all_minus90 = pd.concat([WP05, WP10, WP20, CA05, CA10, CA20, SB05, SB10, SQ12])
    all_plus90.columns = ['X', 'x', 'y', 'z']
    all_minus90.columns = ['X', 'x', 'y', 'z']
    all_plus90_dict = all_plus90.set_index('X').T.to_dict('list')
    all_minus90_dict = all_minus90.set_index('X').T.to_dict('list')

    return all_plus90_dict, all_minus90_dict


def rotate(coords, rotation):
    if rotation == "plus90":
        R = np.matrix([[-0.000000043711, -1.000000000000, 0.000000000000, 0.000000000000],
                       [1.000000000000, -0.000000043711, 0.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif rotation == "minus90":
        R = np.matrix([[-0.000000043711, 1.000000000000, 0.000000000000, 0.000000000000],
                       [-1.000000000000, -0.000000043711, 0.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    else:
        print('choose rotation')
    parallel_coords = {}
    for name in coords:
        coords[name].append(1)
        parallel_coords[name] = np.matmul(R, coords[name])
        parallel_coords[name] = parallel_coords[name].tolist()
        parallel_coords[name] = parallel_coords[name][0]

    return parallel_coords


def main():
    print('Import and complie coords')
    all_plus90_dict, all_minus90_dict = import_coords_compile(PATH)

    parallel_plus90 = rotate(all_plus90_dict, "plus90")
    parallel_minus90 = rotate(all_minus90_dict, "minus90")

    parallel_plus90_df = pd.DataFrame(parallel_plus90).T
    parallel_minus90_df = pd.DataFrame(parallel_minus90).T

    all_annotations = pd.concat([parallel_plus90_df, parallel_minus90_df])
    all_annotations.columns = ['x', 'y', 'z', 'p']
    del all_annotations['p']
    print("Exporting coords to results")
    all_annotations.to_csv('{}/all_annotations_X_DEPTH_parallel_copy.txt'.format(PATH))


main()
