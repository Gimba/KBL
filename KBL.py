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
import operator
import sys
from argparse import RawTextHelpFormatter

import matplotlib.pyplot as plt

import Kullback_Leibler_lib as k


def main(args):
    parser = argparse.ArgumentParser(description=r'Calculate the KBL value for the angle files (.xvg format) found in '
                                                 'the two specified directories.\n\n Files that get generated: \n\n{'
                                                 'output_filename} - .pymol that can be loaded in pymol (@{'
                                                 'output_filename} and colors residues with low kbl-value white and '
                                                 'high kbl-values blue\n'
                                                 'hist.{output_filename}.png - image of the generated histogram\n'
                                                 'hist.{outp_filemane}.dat - data of the histogram',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('dir1', help='reference directory, e.g. ./WT/gromacs_psi_phi_chi/')
    parser.add_argument('dir2', help='to compare to directory, e.g. ./R2215A/gromacs_psi_phi_chi/')
    parser.add_argument('-m', dest='mutations', nargs='?', help='mutations introduced to the structure of the second ' \
                                                                'directory. e.g. ARG2215,ALA2215', default="")
    parser.add_argument('-o', dest='output_filename', nargs='?', help='name of .pymol file that gets written ('
                                                                      'default: kbl_{dir1}_{e1}_{dir2}_{e2}_{'
                                                                      'angles}.pymol')
    parser.add_argument('-a', dest='angles', nargs='?', help='specify angles that should be used to calculate kbl, '
                                                             'e.g. psi,phi,chi1',
                        default="")
    parser.add_argument('-e1', dest='end_dir1', nargs='?', help='lines/number of frames to read from files in '
                                                                'dir1, e.g. 4000')
    parser.add_argument('-e2', dest='end_dir2', nargs='?', help='lines/number of frames to read from files in '
                                                                'dir2, e.g. 20000')
    parser.add_argument('-s', dest='show', action='store_true', help='show histogram', default=False)
    parser.add_argument('-r', dest='resids', nargs='?', help='specify a set of residues that are used in the '
                                                             'calculation', default=False)
    parser.add_argument('-y_range', nargs='?', help='lower and upper bound for y-axis, used to make different kbl '
                                                    'calculations comperable')
    parser.add_argument('-suff', dest='suffix', nargs='?', help='add suffix to file names', default=False)

    args = parser.parse_args()

    # set up input parameters
    dir1 = args.dir1.strip('\'')
    dir2 = args.dir2.strip('\'')

    residues = []
    if args.resids:
        residues = list(map(int, args.resids.split(',')))

    y_range = []
    if args.y_range:
        y_range = args.y_range.split(',')
        y_range = list(map(float, y_range))

    if args.output_filename:
        kbl_filename = args.output_filename

    else:
        dir1_name = dir1.split('/')
        dir2_name = dir2.split('/')

        # if the directories do not differ in their last path-dir us the second last dir for naming
        if dir1_name[-2] == dir2_name[-2]:
            kbl_filename = \
                'kbl_' + dir1_name[-3] + '_' + dir2_name[-3] + '_' + args.angles.replace(',', '_') + '.pymol'
        else:
            kbl_filename = \
                'kbl_' + dir1_name[-2] + '_' + dir2_name[-2] + '_' + args.angles.replace(',', '_') + '.pymol'

    if args.suffix:
        kbl_filename = kbl_filename.split('.')
        kbl_filename = kbl_filename[0] + "_" + args.suffix + "." + kbl_filename[1]

    hist_dat_filename = 'hist.' + kbl_filename.replace('.pymol', '.dat')
    hist_png_filename = 'hist.' + kbl_filename.replace(".pymol", '.png')

    mutations = args.mutations.split(',')

    angles = args.angles.split(',')

    end_dir1 = 999999999
    if args.end_dir1:
        end_dir1 = args.end_dir1
        kbl_filename = kbl_filename.split('_')
        kbl_filename = kbl_filename[0] + "_" + kbl_filename[1] + "(1-" + end_dir1 + ")_" + "_".join(kbl_filename[
                                                                                                    2:])
        end_dir1 = int(end_dir1)

    end_dir2 = 999999999
    if args.end_dir2:
        end_dir2 = args.end_dir2
        kbl_filename = kbl_filename.split('_')
        kbl_filename = kbl_filename[0] + "_" + kbl_filename[1] + "_" + kbl_filename[2] + "(1-" + end_dir2 + ")_" + \
                       "_".join(
                           kbl_filename[3:])
        end_dir2 = int(end_dir2)

    print("pymol output-file:", kbl_filename)

    print('reading data: ', dir1, dir2)
    data = k.read_in_data(dir1, dir2, end1=end_dir1, end2=end_dir2, mutations=mutations, angles=angles,
                          residues=residues)

    print('calculating distributions')
    distributions = k.get_distributions(data)

    print('calculating KBL-values')
    kbl = k.get_kbl(distributions)

    # generate histogram
    all_kbl = []
    all_kbl_names = []
    for i in kbl:
        plt.scatter(float(i[3:]), kbl[i], c="g")
        if len(y_range) > 1:
            plt.ylim(y_range)
        all_kbl.append(kbl[i])
        all_kbl_names.append(i)
    plt.xlabel('residue')
    plt.ylabel('kbl value')

    plt.savefig(hist_png_filename)

    if args.show:
        plt.show()

    # KBL by values
    kbl_sorted = sorted(kbl.items(), key=operator.itemgetter(1), reverse=True)

    # write histogram data to file
    with open(hist_dat_filename, 'w') as o:
        for item in kbl_sorted:
            o.write(item[0] + "," + str(item[1]) + "\n")

    # output max and second maximal value to console
    print('maximum: ', kbl_sorted[0][0], kbl_sorted[0][1])
    print('second : ', kbl_sorted[1][0], kbl_sorted[1][1])
    maximum = kbl_sorted[0][1]

    # generate and write .pymol file
    with open(kbl_filename, 'w') as f:
        f.write('hide lines\nshow cartoon\n')
        for j in kbl:
            if maximum == 0:
                print("comparing the same distributions, exiting...")
                exit()
            val = 1. - (kbl[j] / maximum)
            res = ''.join([i for i in j if i.isdigit()])
            f.write('set_color %sred = [%f, %f, 1]\n' % (j, val, val))
            f.write('color %sred, resid %s\n' % (j, res))


if __name__ == '__main__':
    main(sys.argv)
