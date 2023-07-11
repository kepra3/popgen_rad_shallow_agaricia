#!/usr/anaconda3/bin/python/
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 19/3/22
@description: Take annotation points and export rotated and scaled plot.
"""

# Modules
import argparse
import numpy as np
import open3d as o3d
import json
import pandas as pd

IGNORE_ANNOTATIONS = ['left', 'right']


class Viscore_metadata(object):
    """ Defining viscore metadata as a class object"""

    def __init__(self, subsets_filename, short_name):
        subsets = json.load(open('{}/{}'.format(path, subsets_filename)))
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


def calculate_euler_angle(opposite, adjacent):
    """ Formula for calculating euler angles """
    # Formula for:
    # atan2(opposite, adjacent)
    if adjacent > 0:
        theta = np.arctan(opposite / adjacent) * 180 / np.pi
    elif adjacent < 0 & opposite >= 0:
        theta = (np.arctan(opposite / adjacent) + np.pi) * 180 / np.pi
    elif adjacent < 0 & opposite < 0:
        theta = (np.arctan(opposite / adjacent) - np.pi) * 180 / np.pi
    elif adjacent == 0 & opposite > 0:
        theta = (np.arctan(opposite / adjacent) + np.pi / 2) * 180 / np.pi
    elif adjacent == 0 & opposite < 0:
        theta = (np.arctan(opposite / adjacent) - np.pi / 2) * 180 / np.pi
    else:
        theta = None
        print("theta is undefined")
    return theta.__float__()


def rotate_matrix(pcd, up_vector):
    """ Create rotation matrix from euler angles and up vector """
    origin = [0, 0, 0]
    # angle xz, theta
    x_diff = origin[0] - up_vector[0]  # opposite - reef perpendicular
    z_diff = origin[2] - up_vector[2]  # adjacent - depth
    theta_xz = calculate_euler_angle(x_diff, z_diff)
    print('Theta is ...', theta_xz)
    # angle yz, psi
    y_diff = origin[1] - up_vector[1]  # opposite - reef parallel, y-coord stays the same
    z_diff = origin[2] - up_vector[2]  # adjacent - depth, z-coord is changed
    psi_yz = calculate_euler_angle(y_diff, z_diff)
    print('Psi is ...', psi_yz)
    # needs radians input
    theta_xz_radians = theta_xz / 180 * np.pi
    psi_yz_radians = psi_yz / 180 * np.pi
    R = pcd.get_rotation_matrix_from_xyz((psi_yz_radians, theta_xz_radians, 0))
    return R


def main(annotations_filename, subsets_filename, short_name, path):
    """ Rotating annotations """
    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)
    print('Read assignment file ...')
    annotations = get_annotations('{}/{}'.format(path, annotations_filename))
    print('Create rotation matrix ...')
    up_vector = viscore_md.dd[0:3]
    # Make up a point cloud in order to access pcd rotation formulae
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector([[1, 1, 1]])
    R = rotate_matrix(pcd, up_vector)
    print('Rotate and scale annotations ...')
    rotated_annotations = {}
    for name in annotations:
        rotated_annotations[name] = np.matmul(R, annotations[name])
        rotated_scaled_annotations = rotated_annotations
        for i in range(3):
            rotated_scaled_annotations[name][i] = rotated_annotations[name][i] * viscore_md.scale_factor

    print('Exporting rotated and scaled annotations as csv')
    annotation_df = pd.DataFrame(rotated_scaled_annotations).T
    annotation_df.columns = ['x', 'y', 'z']
    annotation_df.to_csv('../results/rotated_scaled_annotations_X_{}_copy.csv'.format(short_name))

if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Rotate annotations X")
    parser.add_argument('ply_filename', help='Filename of PLY file')
    parser.add_argument('annotations_filename', help='Filename of TXT file')
    args = parser.parse_args()

    ply_filename = args.ply_filename + ".ply"
    annotations_filename = args.annotations_filename + ".txt"
    short_name = "_".join(ply_filename.split('_')[0:4])
    subsets_filename = short_name + "_subsets.json"
    path = "../data"

    # e.g.,
    # args.ply_filename = "cur_kal_20m_20200214_decvis_02"

    main(annotations_filename, subsets_filename, short_name, path)
