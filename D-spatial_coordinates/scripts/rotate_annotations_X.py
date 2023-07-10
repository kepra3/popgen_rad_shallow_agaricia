#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 19/3/22
@description: Take plot point and annotation points and export rotated and scaled plot
 and choice of whether building a KDTree.
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


def fit_a_plane_ransac(pcd):
    """ Fit a plane to a point cloud using the RANSAC method """
    plane_model, inliers = pcd.segment_plane(distance_threshold=0.001,
                                             ransac_n=3,
                                             num_iterations=1000)
    [a, b, c, d] = plane_model
    print(f"Plane equation: {a:.2f}x + {b:.2f}y + {c:.2f}z + {d:.2f} = 0")
    return plane_model, inliers


def calc_plane_angles(plane_model):
    """ Calculate the angles of a plane through the plane equation
        Angles are calculated by finding the angle between two vectors """
    plane_normal = [plane_model[0], plane_model[1], plane_model[2]]
    slope_xz = plane_normal[0] / plane_normal[2]
    slope_yz = plane_normal[1] / plane_normal[2]
    theta = np.arctan(slope_xz) * 180 / np.pi
    print('The angle between x and z is ...', theta)  # about y-axis (xz)
    psi = np.arctan(slope_yz) * 180 / np.pi
    print('The angle between y and z is ...', psi)  # about x-axis (yz)
    # the angle between the x-y plane... i.e., the elevation
    xy_normal = [0, 0, 1]
    mag_xy = np.linalg.norm(xy_normal)
    mag_plane = np.linalg.norm(plane_normal)
    cos_elevation = np.dot(xy_normal, plane_normal) / (mag_xy * mag_plane)
    elevation = np.arccos(cos_elevation) * 180 / np.pi
    print('the angle between plane and x-y plane is ...', elevation)
    return theta.__float__(), psi.__float__(), elevation.__float__()


def main(ply_filename, annotations_filename, subsets_filename, path, KDTree='No'):
    """ Write description here """

    short_name = "_".join(ply_filename.split('_')[0:4])
    # Create plot result file
    plot_info = '../results/plot_info.txt'
    with open(plot_info, 'a') as results_out:
        if results_out.tell() == 0:
            print("Creating a new file\n")
            results_out.write("plot_name\tplot_points\ttheta_xz\tpsi_yz\televation\n")
        else:
            "File exists"

    # 1. PREPARATION SUBSET COLONY POINTS AND SCALE & ROTATE ALL POINTS ####
    print('Reading PLY file {} ...'.format(ply_filename))
    pcd = o3d.io.read_point_cloud('{}/{}'.format(path, ply_filename))
    print('Read viscore metadata file ...')
    viscore_md = Viscore_metadata(subsets_filename, short_name)
    print('Rotating matrix ...')
    up_vector = viscore_md.dd[0:3]
    R = rotate_matrix(pcd, up_vector)
    #pcd_r = copy.deepcopy(pcd)
    #pcd_r.rotate(R, center=(0, 0, 0))
    #print('Scaling point cloud ...')
    #pcd_r.scale(viscore_md.scale_factor, center=(0, 0, 0))
    print('Read assignment file ...')
    annotations = get_annotations('{}/{}'.format(path, annotations_filename))
    print('Rotate and scale annotations ...')
    rotated_annotations = {}
    for name in annotations:
        rotated_annotations[name] = np.matmul(R, annotations[name])
        rotated_scaled_annotations = rotated_annotations
        for i in range(3):
            rotated_scaled_annotations[name][i] = rotated_annotations[name][i] * viscore_md.scale_factor
    #print('Get scaled ranges for each sample...')
    #ranges = get_ranges('{}/{}'.format(path, annotations_filename), annotations, viscore_md.scale_factor)


    #Write info about plot
    #print('Get angles of construct...')
    #plane_model, inliers = fit_a_plane_ransac(pcd_r)
    #theta_xz, psi_yz, elevation = calc_plane_angles(plane_model)
    #print('Write plot angles to file')
    #with open(plot_info, 'a') as results_out:
    #    results_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(short_name, len(np.asarray(pcd_r.points)),
    #                                                         theta_xz, psi_yz, elevation))
    print('Exporting rotated and scaled annotations as csv')
    annotation_df = pd.DataFrame(rotated_scaled_annotations).T
    annotation_df.columns = ['x', 'y', 'z']
    #annotation_df['range'] = ranges.values()
    annotation_df.to_csv('../results/rotated_scaled_annotations_X_{}.csv'.format(short_name))
    #print('Export ply file to HardDrive')
    #o3d.io.write_point_cloud('/Volumes/KP3/coralscape/{}_rotated_scaled.ply'.format(short_name), pcd_r)
    #if KDTree == 'Yes':
    #    print('Building KDTree ...')
    #    pcd_tree_r = o3d.geometry.KDTreeFlann(pcd_r)
    #    # Save KDTree
    #    with open('/Volumes/KP3/coralscape/{}_KDTree.pkl'.format(short_name), "wb") as kdtree_file:
    #        pickle.dump(pcd_tree_r, kdtree_file)
    #else:
    #    'Not building KDTree'


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="Rotate annotations X")
    parser.add_argument('ply_filename', help='Filename of PLY file')
    parser.add_argument('annotations_filename', help='Filename of annotations file')
    parser.add_argument('subsets_filename', help='Filename of Viscore metadata, subsets.json')
    parser.add_argument('site', help='Site code for location')
    parser.add_argument('KDTree', help='Whether you want to build and export a KDTree or not')
    args = parser.parse_args()

    ply_filename = args.ply_filename
    annotations_filename = args.annotations_filename
    subsets_filename = args.subsets_filename

    path = "/Users/kprata/Dropbox/agaricia_project_2019/shalo_ag/Photogrammetry/CloudCompare/{}".format(args.site)
    KDTree = args.KDTree
    # e.g.,
    # args.ply_filename = "cur_kal_20m_20200214_decvis_02.ply"
    # args.annotations_filename = "cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt"
    # args.subsets_filename = "cur_kal_20m_20200214_subsets.json"
    # args.site = "WP20"
    # args.KDTree = "No"

    # PLOTS & ROTATIONS:
    # WP05 cur_kal_05m_20200214_decvis_02_KP.txt
    # theta -6.91, psi 21.08
    # WP10 cur_kal_10m_20200214_decvis_02_KP_905_updated16-3-22.txt
    # theta = -25.11, psi = 11.65
    # WP20 cur_kal_20m_20200214_decvis_02_KP_16-12-21_completed.txt
    # theta = 9.02, psi = 19.25
    # SB05 cur_sna_05m_20200303_decvis_02_SF_HI_19-1-22.txt
    # theta = -3.90, psi = 12.30
    # SB10 cur_sna_10m_20200303_decvis_02
    # SB10 cur_sna_10m_20201202_decvisann_HI_14-12-21.txt
    # theta = -0.72, psi = -2.71
    # SB20 cur_sna_20m_20200303_decvis_02
    # SB20 cur_sna_20m_20190410_decvisann_HI_12_12.txt
    # theta = -16.82, psi = 18.23
    # CA05 cur_cas_05m_20201212_decvis_02_KP_31-01-22.txt
    # note had to alter subsets.json file to have cur_cas_05m_20201212/cur_cas_05m_20201212
    # theta = 2.86, psi = -0.24
    # CA10 cur_cas_10m_20201212_decvis_02.ply cur_cas_10m_20201210_decvis_02_SH_done.txt
    # note had to alter subsets.json
    # Theta is ... -12.36016759332716
    # Psi is ... -10.45657839852049
    # CA20 cur_cas_20m_20201212_decvis_02.ply cur_cas_20m_20201212_decvis_02_SH_01-06-22.txt
    # NOTE 1161 did not have a left (thus have put an _X for now - need to fix this)
    # Theta is ... -32.03105065955869
    # Psi is ... -6.77682359250437
    # SQ12 cur_seb_10m_20201210_decvis_02.ply cur_seb_10m_20201210_decvis_02_SH_updated.txt
    #Theta is ... -11.644224255407794
    #Psi is ... 2.8498421712287363
    # SQ20 cur_seb_20m_20201210_decvis_02.ply cur_seb_20m_20201210_decvis_02_SH_10-02-2022.txt
    # Theta is ... 29.32225489795626
    # Psi is ... -0.6445221821626353

    main(ply_filename, annotations_filename, subsets_filename, path, KDTree)
