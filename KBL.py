#! /usr/bin/env python

# Copyright (c) 2019 Martin Rosellen

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import math
import operator
import sys
from argparse import RawTextHelpFormatter
import os
import matplotlib
from re import findall

matplotlib.use('agg')

import matplotlib.pyplot as plt

import Kullback_Leibler_lib as k
import time

trajectory_filetypes = [".nc"]


def main(args):
    tic = time.time()
    parser = argparse.ArgumentParser(
        description='Calculate the Jensen Shannon divergence between angle distributions extracted from specified '
                    'files\n\n'
                    'Files that get generated \n\n{'
                    'output_filename} - .pml that can be loaded in pml (@{'
                    'output_filename} and colors residues with low kbl-value white and '
                    'high kbl-values blue\n'
                    'hist.{output_filename}.png - image of the generated histogram\n'
                    'hist.{outp_filemane}.dat - data of the histogram',
        formatter_class=RawTextHelpFormatter)
    # TODO: give option to read in whole directories and specify how many frames should be considered
    parser.add_argument('files1', help='list of trajectory files (comma separated) or directory')
    parser.add_argument('files2', help='list of trajectory files (comma separated) or directory')
    parser.add_argument('-m', dest='mutations', nargs='?', help='mutations introduced to the structure of the second ' \
                                                                'directory. e.g. ARG2215,ALA2215', default="")
    parser.add_argument('-o', dest='output_filename', nargs='?', help='name of .pml file that gets written ('
                                                                      'default: kbl_{dir1}_{e1}_{dir2}_{e2}_{'
                                                                      'angles}.pml')
    parser.add_argument('-a', dest='angles', nargs='?', help='specify angles that should be used to calculate kbl, '
                                                             'e.g. psi,phi,chi1',
                        default="")
    parser.add_argument('-e1', dest='end_dir1', nargs='?', help='lines/number of frames to read from files in '
                                                                'dir1, e.g. 4000')
    parser.add_argument('-e2', dest='end_dir2', nargs='?', help='lines/number of frames to read from files in '
                                                                'dir2, e.g. 20000')
    parser.add_argument('-s', dest='show', action='store_true', help='show histogram', default=False)
    parser.add_argument('-r', dest='resids', nargs='?', help='specify a set of residues that are used in the '
                                                             'calculation, default: all residues', default=False)
    parser.add_argument('-y_range', nargs='?', help='lower and upper bound for y-axis, used to make different kbl '
                                                    'calculations comperable')
    parser.add_argument('-suff', dest='suffix', nargs='?', help='add suffix to file names', default=False)
    parser.add_argument('-t1', dest='topology1', nargs='?', help='topology for files1')
    parser.add_argument('-t2', dest='topology2', nargs='?', help='topology for files2')
    args = parser.parse_args()

    # set up input parameters
    files1 = args.files1
    files2 = args.files2

    # files given as directory
    if files1[-1] == "/":
        print("Directory with trajectories (files1): ", files1)
        files1 = [files1 + f for f in os.listdir(files1) if
                  f.split(".")[-1] in trajectory_filetypes[0] and "prod_" in f]
        # sort files, index expected to be last number in file name
        files1.sort(key=lambda x: int(findall(r'\d+', x)[-1]))
        print("Found files: ", files1)

    # files given as a list of files
    elif "," in files1:
        files1 = files1.split(",")
        print("Trajectory files (files1): ", files1)
    # single file
    else:
        files1 = [files1]
        print("Trajectory file (files1): ", files1[0])

    # files given as directory
    if files2[-1] == "/":
        print("Directory with trajectories (files2): ", files2)
        files2 = [files2 + f for f in os.listdir(files2) if
                  f.split(".")[-1] in trajectory_filetypes[0] and "prod_" in f]
        files2.sort(key=lambda x: int(findall(r'\d+', x)[-1]))
        print("Found files: ", files2)

    # files given as a list of files
    elif "," in files2:
        files2 = files2.split(",")
        print("Trajectory files (files2): ", files2)
    # single file
    else:
        files2 = [files2]
        print("Trajectory file (files2): ", files2[0])

    residues = []
    if args.resids:
        if "," in args.resids:
            residues = list(map(int, args.resids.split(',')))
        elif "-" in args.resids:
            start, end = args.resids.split("-")
            residues = [l for l in range(int(start), int(end) + 1)]
        else:
            residues = [args.resids]
            residues = list(map(int, residues))

    y_range = []
    if args.y_range:
        y_range = args.y_range.split(',')
        y_range = list(map(float, y_range))

    if args.output_filename:
        kbl_filename = args.output_filename

    else:
        dir1_name = files1.split('/')
        dir2_name = files2.split('/')

        # if the directories do not differ in their last path-dir us the second last dir for naming
        # e.g.: somedir/naming_dir1/subdir; somedir/naming_dir2/subdir
        if dir1_name[-2] == dir2_name[-2]:
            kbl_filename = \
                'kbl_' + dir1_name[-3] + '_' + dir2_name[-3] + '_' + args.angles.replace(',', '_') + '.pml'
        else:
            kbl_filename = \
                'kbl_' + dir1_name[-2] + '_' + dir2_name[-2] + '_' + args.angles.replace(',', '_') + '.pml'

    if args.suffix:
        kbl_filename = kbl_filename.split('.')
        kbl_filename = kbl_filename[0] + "_" + args.suffix + "." + kbl_filename[1]

    # filenames for histogram and its data (data sorted descending by kbl value)
    hist_dat_filename = 'hist.' + kbl_filename.replace('.pml', '.dat')
    hist_png_filename = 'hist.' + kbl_filename.replace(".pml", '.png')

    mutations = args.mutations.split(',')

    angles = args.angles.split(',')
    n_frames1 = 9999999999
    n_frames2 = 9999999999
    if args.end_dir1:
        n_frames1 = int(args.end_dir1)
    if args.end_dir2:
        n_frames2 = int(args.end_dir2)

    print("pymol output-file:", kbl_filename)
    print('topologies:',args.topology1, args.topology2)
    print('reading data:', files1, files2)
    # data = k.read_in_data(files1, files2, end1=end_dir1, end2=end_dir2, mutations=mutations, angles=angles,
    #                       residues=residues, topologies=[args.topology1, args.topology2])
    data = k.read_in_data(files1, files2, mutations=mutations, angles=angles, residues=residues,
                          topologies=[args.topology1, args.topology2], n_frames1=n_frames1, n_frames2=n_frames2)
    print('calculating distributions')
    distributions = k.get_distributions(data)

    print('calculating JSD-values')
    jsd = k.get_jsd(distributions)

    # generate histogram

    if args.show:
        # TODO: give option to plot only a range of residues (the 10 highest JSD values)
        all_jsd = []
        all_jsd_names = []
        for name, distance in jsd.items():
            plt.scatter(name, distance, c="g")
            if len(y_range) > 1:
                plt.ylim(y_range)
            all_jsd.append(distance)
            all_jsd_names.append(name)
        plt.xlabel('residue')
        plt.ylabel('kbl value')

        plt.savefig(hist_png_filename)
        plt.show()

    # sort KBL values for saving in hist data file
    flatter_jsd_dict = {}
    for l, d in jsd.items():
        for u, f in d.items():
            flatter_jsd_dict[l + " " + u] = f

    kbl_sorted = sorted(flatter_jsd_dict.items(), key=operator.itemgetter(1), reverse=True)

    # write histogram data to file
    with open(hist_dat_filename, 'w') as o:
        for item in kbl_sorted:
            o.write(item[0] + "," + str(item[1]) + "\n")

    # output max and second maximal value to console
    print('maximum: ', kbl_sorted[0][0], kbl_sorted[0][1])
    print('second : ', kbl_sorted[1][0], kbl_sorted[1][1])
    # maximum = kbl_sorted[0][1]

    maxima = sorted(flatter_jsd_dict.items(), key=operator.itemgetter(0), reverse=True)

    maxis = []
    for l in range(0, len(maxima) - 1, 2):
        maxis.append((maxima[l][1] + maxima[l + 1][1]))
    psi_phi_max = max(maxis)

    color_jsd_dict = {}

    # for key1 in list(jsd.keys()):
    #     for key2 in list(jsd.keys()):
    #         if key1.split(" ")[0] in key2:
    #             if not key1.split(" ")[0] in color_jsd_dict.keys():
    #                 color_jsd_dict[key1.split(" ")[0]] = {}
    #             temp_dict = {key2.split(" ")[1]: jsd[key2]}
    #             color_jsd_dict[key1.split(" ")[0]].update(temp_dict)
    #             del jsd[key2]
    # generate and write .pml file
    colors = []
    with open(kbl_filename, 'w') as f:
        f.write('hide lines\nshow cartoon\n')
        f.write('# ' + str(kbl_sorted))
        for resid, angles in jsd.items():
            if "phi" in angles.keys():
                phi_val = angles["phi"]
            else:
                phi_val = 0
            if "psi" in angles.keys():
                psi_val = angles["psi"]
            else:
                psi_val = 0

            backbone = (psi_val + phi_val)/ 2
            colors.append(backbone)

            color_jsd_dict[resid] = backbone

            # elif backbone < 0.5:
            #     backbone = 0
            if "chi" in angles.keys():
                chi_val = math.sqrt(angles["chi"])
            else:
                chi_val = 0
        colors.sort()
        min_colors = min(colors)
        max_colors = max(colors)
        max2_colors = colors[-2]
        for resid, color in color_jsd_dict.items():
            res = ''.join([i for i in resid if i.isdigit()])
            if color == max_colors:
                f.write('label resid ' + res + ' and n. CB,  \'' + resid + '=' + str(round(color, 2)) + '\'\n')
            if color == max2_colors:
                f.write('label resid ' + res + ' and n. CB,  \'' + resid + '=' + str(round(color, 2)) + '\'\n')
            color = (color - min_colors) / (max_colors - min_colors)
            f.write('set_color %sred = [%f, %f, %f]\n' % (resid, 1, 1 - color, 1 - color))
            f.write('color %sred, resid %s\n' % (resid, res))
    toc = time.time()
    print("------ ", round(toc - tic, 4), "seconds ------")


if __name__ == '__main__':
    main(sys.argv)
