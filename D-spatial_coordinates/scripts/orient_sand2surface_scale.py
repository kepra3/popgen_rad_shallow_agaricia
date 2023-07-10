#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 3/5/22
@description: Orient plot to the XY-axis plane.
"""

import argparse
import json
# import open3d as o3d
import numpy as np
# import copy
import pandas as pd


class Viscore_metadata(object):
    """ Defining viscore metadata as a class object"""

    def __init__(self, subsets_filename, short_name):
        subsets = json.load(open('{}/{}'.format(PATH, subsets_filename)))
        # Determine json key of primary model
        if '{0}/{0}'.format(short_name) in subsets['d']:
            subsets_ortho = subsets['d']['{0}/{0}'.format(short_name)]['c']['ortho']
        elif self.short_name in subsets['d']:
            subsets_ortho = subsets['d'][short_name]['c']['ortho']
        else:
            print('Model not found in subsets.json!')
        self.dd = subsets_ortho['dd']
        self.scale_factor = subsets_ortho['scale_factor']
        self.r = subsets_ortho['vecs']['r']
        self.u = subsets_ortho['vecs']['u']
        self.n = subsets_ortho['vecs']['n']
        self.c = subsets_ortho['vecs']['c']
        self.cc = subsets_ortho['vecs']['cc']
        self.cam_up = subsets_ortho['vecs']['cam']['up']
        self.cam_eye = subsets_ortho['vecs']['cam']['eye']
        self.cam_target = subsets_ortho['vecs']['cam']['target']


def get_annotations(annotations_path):
    """ Read annotations from txt file """
    annotations = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        if any(flag in line for flag in IGNORE_ANNOTATIONS):
            continue
        cols = line.rstrip().replace(',', ' ').split()
        annotations[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    return annotations


def get_ranges(annotations_path, annotations, scale):
    """ Get longest length of colony divided by 2 """
    complete_set = {}
    annotations_file = open(annotations_path, 'r')
    for line in annotations_file:
        cols = line.rstrip().replace(',', ' ').split()
        complete_set[cols[0]] = [float(i) for i in cols[1:4]]
    annotations_file.close()
    ranges = {}
    for name in annotations:
        # euclidean distance
        x1 = complete_set['{}_left'.format(name)][0]
        x2 = complete_set['{}_right'.format(name)][0]
        y1 = complete_set['{}_left'.format(name)][1]
        y2 = complete_set['{}_right'.format(name)][1]
        z1 = complete_set['{}_left'.format(name)][2]
        z2 = complete_set['{}_right'.format(name)][2]
        ranges[name] = (((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5 / 2) * scale
    return ranges


def main(ply_filename, annotations_filename, PATH):
    """ Orienting substrate parallel to surface and scale xyz """
    short_name = "_".join(ply_filename.split('_')[0:4])
    subsets_filename = "_".join([short_name, 'subsets.json'])
    # print('Reading PLY file {} ...'.format(ply_filename))
    # pcd = o3d.io.read_point_cloud('{}/{}.ply'.format(PATH, ply_filename))
    print('Read viscore metadata file ...')
    if short_name == "cur_cas_10m_20201212":
        scale = 0.10383552972538619
    elif short_name == "cur_cas_20m_20201212":
        scale = 0.10214017812800297
    elif short_name == "cur_seb_10m_20201210":
        scale = 0.21249641478078696
    elif short_name == "cur_seb_20m_20201210":
        scale = 0.1937345919527237
    else:
        viscore_md = Viscore_metadata(subsets_filename, short_name)
        scale = viscore_md.scale_factor
    print('Read annotations file ...')
    annotations = get_annotations('{}/{}'.format(PATH, annotations_filename))
    if short_name == 'cur_kal_05m_20200214':
        R = np.matrix([[0.999081313610, -0.007136042230, -0.042256940156, -0.222624540329],
                       [-0.007136042230, 0.944570958614, -0.328229755163, -1.729230880737],
                       [0.042256940156, 0.328229755163, 0.943652272224, 5.980512142181],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_kal_10m_20200214':
        R = np.matrix([[0.972240805626, -0.043610032648, -0.229882493615, -8.932192802429],
                       [-0.043610032648, 0.931488096714, -0.361148327589, -14.032592773438],
                       [0.229882493615, 0.361148327589, 0.903728902340, -34.574745178223],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_kal_20m_20200214':
        R = np.matrix([[0.997460484505, -0.011900402606, -0.070220448077, -0.287017107010],
                       [-0.011900402606, 0.944233059883, -0.329062759876, -1.345001220703],
                       [0.070220448077, 0.329062759876, 0.941693544388, 9.885766983032],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_cas_05m_20201212':
        R = np.matrix([[0.994036436081, -0.001632568543, -0.109036169946, -0.633613109589],
                       [-0.001632568543, 0.999553084373, -0.029849467799, -0.173454284668],
                       [0.109036169946, 0.029849467799, 0.993589520454, 1.386008739471],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_cas_10m_20201212':
        R = np.matrix([[0.978240072727, 0.024879366159, -0.205978974700, -4.229547023773],
                       [0.024879366159, 0.971553981304, 0.235507741570, 4.835884094238],
                       [0.205978974700, -0.235507741570, 0.949794054031, -10.681835174561],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_cas_20m_20201212':
        R = np.matrix([[0.999715030193, 0.000036135007, -0.023871334270, -0.431133270264],
                       [0.000036135007, 0.999995410442, 0.003027042607, 0.054670333862],
                       [0.023871334270, -0.003027042607, 0.999710440636, 0.221694946289],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_sna_05m_20200303':
        R = np.matrix([[0.999736130238, -0.001426698873, -0.022926803678, -0.151225090027],
                       [-0.001426698873, 0.992286145687, -0.123960413039, -0.817646026611],
                       [0.022926803678, 0.123960413039, 0.992022275925, -2.414521217346],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_sna_10m_20200303':
        R = np.matrix([[0.949437201023, 0.001126733376, -0.313954859972, -2.941060304642],
                       [0.001126733376, 0.999974846840, 0.006996125914, 0.065540313721],
                       [0.313954859972, -0.006996125914, 0.949412107468, -0.251988410950],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_sna_20m_20200303':
        R = np.matrix([[0.950774013996, 0.045413658023, 0.306539326906, -0.731462478638],
                       [0.045413658023, 0.958103418350, -0.282799273729, 0.674812316895],
                       [-0.306539326906, 0.282799273729, 0.908877432346, 24.230812072754],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_seb_10m_20201210':
        R = np.matrix([[0.952747523785, -0.024975303560, -0.302734822035, -1.868291378021],
                       [-0.024975303560, 0.986799299717, -0.160010576248, -0.987485885620],
                       [0.302734822035, 0.160010576248, 0.939546823502, 0.715722084045],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    elif short_name == 'cur_seb_20m_20201210':
        R = np.matrix([[0.947396755219, -0.008276428096, 0.319954514503, 2.978732347488],
                       [-0.008276428096, 0.998697817326, 0.050340626389, 0.468662261963],
                       [-0.319954514503, -0.050340626389, 0.946094572544, 0.607359886169],
                       [0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000]])
    else:
        print("Issue with filename")
    # print('Rotating ...')
    # pcd_r = copy.deepcopy(pcd)
    # pcd_r.rotate(R, center=(0, 0, 0))
    # print('Scaling point cloud ...')
    # pcd_r.scale(viscore_md.scale_factor, center=(0, 0, 0))
    print('Rotate and scale annotations ...')
    rotated_annotations = {}
    for name in annotations:
        annotations[name].append(1)
        rotated_annotations[name] = np.matmul(R, annotations[name])
        scaled_rotated_annotations = rotated_annotations
        scaled_rotated_annotations[name] = scaled_rotated_annotations[name].tolist()
        scaled_rotated_annotations[name] = scaled_rotated_annotations[name][0]
        for i in range(3):
            scaled_rotated_annotations[name][i] = rotated_annotations[name][i] * scale
    print("Scaling ranges for each sample")
    #ranges = get_ranges('{}/{}'.format(PATH, annotations_filename), annotations, scale)
    print('Exporting rotated and scaled annotations as csv')
    annotation_df = pd.DataFrame(scaled_rotated_annotations).T
    annotation_df.columns = ['x', 'y', 'z', 'p']
    # annotation_df['range'] = ranges.values()
    # TODO: NOTE! not including range for now
    annotation_df.to_csv('{}/scaled_HORIZ_annotations_X_{}.csv'.format(PATH, short_name))


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Colony clean and measure")
    parser.add_argument('ply_filename', help='Filename of PLY file')
    parser.add_argument('annotations_filename', help='Filename of annotations file')
    args = parser.parse_args()

    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename

    name_split = ply_filename.split("_")

    if name_split[1] == "kal":
        loc = "WP"
    elif name_split[1] == "cas":
        loc = "CA"
    elif name_split[1] == "sna":
        loc = "SB"
    elif name_split[1] == "seb":
        loc = "SQ"
    else:
        print("issue with ply filename")

    if name_split[2] == "05m":
        depth = "05"
    elif name_split[2] == "10m":
        depth = "10"
    elif name_split[2] == "20m":
        depth = "20"
    else:
        print("issue with ply filename")
    loc_depth = [loc, depth]
    site = "".join(loc_depth)

    # GLOBAL VARIABLES
    IGNORE_ANNOTATIONS = ['left', 'right']
    PATH = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/{}".format(site)

    main(ply_filename, annotations_filename, PATH)
    # PLOTS & ROTATIONS:
    # WP05
    # python orient_sand2surface_scale.py cur_kal_05m_20200214_decvis_02.ply cur_kal_05m_20200214_decvis_02_KP.txt
    # theta -6.91, psi 21.08
    # WP10
    # python orient_sand2surface_scale.py cur_kal_10m_20200214_decvis_02.ply cur_kal_10m_20200214_decvis_02_KP_905_updated16-3-22.txt
    # theta = -25.11, psi = 11.65
    # WP20
    # python orient_sand2surface_scale.py cur_kal_20m_20200214_decvis_02.ply cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt
    # theta = 9.02, psi = 19.25
    # SB05
    # python orient_sand2surface_scale.py cur_sna_05m_20200303_decvis_02.ply cur_sna_05m_20200303_decvis_02_SF_HI_19-1-22.txt
    # theta = -3.90, psi = 12.30
    # SB10
    # python orient_sand2surface_scale.py cur_sna_10m_20200303_decvis_02.ply cur_sna_10m_20201202_decvisann_HI_14-12-21.txt
    # theta = -0.72, psi = -2.71
    # SB20
    # python orient_sand2surface_scale.py cur_sna_20m_20200303_decvis_02.ply cur_sna_20m_20190410_decvisann_HI_12_12.txt
    # theta = -16.82, psi = 18.23
    # CA05
    # python orient_sand2surface_scale.py cur_cas_05m_20201212_decvis_02.ply cur_cas_05m_20201212_decvis_02_KP_31-01-22.txt
    # note had to alter subsets.json file to have cur_cas_05m_20201212/cur_cas_05m_20201212
    # theta = 2.86, psi = -0.24
    # CA10
    # python orient_sand2surface_scale.py cur_cas_10m_20201212_decvis_02.ply cur_cas_10m_20201210_decvis_02_SH_done.txt
    # CA20
    # python orient_sand2surface_scale.py cur_cas_20m_20201212_decvis_02.ply cur_cas_20m_20201212_decvis_02_SH_26-05-22.txt
    # SQ12
    # python orient_sand2surface_scale.py cur_seb_10m_20201210_decvis_02.ply cur_seb_10m_20201210_decvis_02_SH_updated.txt
    # SQ20
    # python orient_sand2surface_scale.py cur_seb_20m_20201210_decvis_02.ply cur_seb_20m_20201210_decvis_02_SH_10-02-2022.txt

    # TODO: check, 0 ranges from _left and _right WP05
    # TODO: KP0160_AC_SQ12 does not have a two longest edge points! Need to find other point
    # TODO: KP1161_AC_CA20 does not have _left
    # TODO: fix script so that when _X, the range = 0