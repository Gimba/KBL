# Copyright (c) 2019 Martin Rosellen, Lucas Siemons

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

import numpy as np
import glob
import pytraj as pt

def read_1_set(path: "", end: "", res: "", angles: list, excluded_angles: list = ['chi1'], excluded_types: list = [
    'GLY', 'ALA']) -> list:
    """
    
    :param path: base path of files
    :param end: read until this line
    :param res: residue to read in e.g. PRO123
    :param angles: angles to include e.g. ['phi','psi','chi1]
    :param excluded_angles: angles that need special attention e.g. chi1 not present for ALA and GLY
    :param excluded_types: do not read in data for angles that have this type
    :return: list with angle values
    """

    res_type = ''.join([i for i in res if not i.isdigit()])

    # generate list with file names
    trajectories = ['productions/prod_1.nc']

    # for angle in angles:
    #     ## leave out files for specified angles and residues
    #     if not (angle in excluded_angles and res_type in excluded_types):
    #         files.append(path + angle + res + '.xvg')
    #     else:
    #         print('omitting file', path + angle + res + '.xvg')
    # # read in data from files
    file_tuple = []
    #     f = open(file)
    #     file_tuple.append(f.readlines())
    #     f.close()
    # file_tuple = tuple(file_tuple)

    # extract angle values from file contents
    all_angles = []
    for i in zip(*file_tuple):
        if i[0][0] not in ('#', '@'):
            if float(i[0].split()[0]) < end:
                angles = [a.split()[1] for a in i]
                angles = [float(a) for a in angles]
                all_angles.append(angles)

    return all_angles

def extract_angles(files, residues, angles):
    """

    :param files: trajectories
    :param residues: range of residues given in cpptraj convention
    :param angles:
    :return:
    """
    # "-1" needed for conversion between python convention (starts with 0) and cpptraj convention (starts with 1)
    residues = [int(r[3:]) - 1 for r in residues]
    traj = pt.TrajectoryIterator(files, "native_hex.prmtop")
    data = pt.multidihedral(traj, dihedral_types=' '.join(angles), resrange=residues)
    # convert data to dictionnary
    data_dict = {}
    for k in data.keys():
        # TODO: map keys back to residue names
        data_dict[k] = list(data[k].values)
    return(data_dict)

def get_residue_names_from_file(dir: "", angle: list) -> list:
    traj = pt.load('native_hex.inpcrd', 'native_hex.prmtop')
    traj.strip(":WAT,CL-")
    resids = []
    for t in traj.top.residues:
        resids.append(t.name + str(t.original_resid))
    return resids


def get_resids_from_files(dir1: "", dir2: "", angles: list, residues: list = []) -> list:
    """
    Returns a list of residue ids extracted from file names present in both directories. Also a log file with 
    residue ids missing in one or the other directory gets written to resids.log 
    
    :param dir1: directory of reference angle values
    :param dir2: directory of the other set of angles
    :param angles: angles consider
    :return: list of residue ids present in both directories
    """

    # create dictionnaries to check whether both directories contain the same residue ids for a specific angle
    resids_dir1 = {}
    resids_dir2 = {}
    angle_resid_intersects = {}

    for angle in angles:
        resids_dir1[angle] = get_residue_names_from_file(dir1, angle)
        resids_dir2[angle] = get_residue_names_from_file(dir2, angle)
        angle_resid_intersects[angle] = set(resids_dir1[angle]) & set(resids_dir2[angle])

        if residues:
            temp_lst = []
            for res in list(angle_resid_intersects[angle]):
                if int(res[3:]) in residues:
                    temp_lst.append(res)
            angle_resid_intersects[angle] = set(temp_lst)

    return angle_resid_intersects


def read_in_data(dir1: "", dir2: "", end1: int = 999999999, end2: int = 999999999, mutations: list = [],
                 angles: list = [], residues: list = []) -> dict:
    angle_mutual_residues = get_resids_from_files(dir1, dir2, angles, residues)

    data = {}

    def add_mutation_res():
        dir1_angles = np.array(read_1_set(dir1, end1, mutations[0], angles)).T
        dir2_angles = np.array(read_1_set(dir2, end2, mutations[1], angles)).T

        # handle case where the mutation has not got the chi1 angle
        if len(dir1_angles) > len(dir2_angles):
            dir1_angles = np.delete(dir1_angles, -1, 0)
        # handle case where the non-mutation has not got the chi1 angle
        elif len(dir1_angles) < len(dir2_angles):
            dir2_angles = np.delete(dir2_angles, -1, 0)
        # naming according to reference residue ids (dir1 pdb)
        data[mutations[0]] = (dir1_angles, dir2_angles)

    # # add specified mutation residue which has different file names
    if len(mutations) > 1:
        if residues:
            # check if given mutation is in specified residues
            if int(mutations[0][3:]) in residues:
                add_mutation_res()

        else:
            add_mutation_res()

    # for res in angle_mutual_residues[angles[0]]:
    dir1_angles = extract_angles(["productions/prod_1.nc"], angle_mutual_residues[angles[0]], angles)

    dir2_angles = extract_angles(["productions/prod_1.nc"], angle_mutual_residues[angles[0]], angles)
    # TODO: adjust object types to match the ones used before pytraj
    # dir1_angles = np.array(read_1_set(dir1, end1, res, angles)).T
    # dir2_angles = np.array(read_1_set(dir2, end2, res, angles)).T
    data = (dir1_angles, dir2_angles)

    return data


def make_hist(data):
    not_chi1 = ('GLY', 'ALA')

    bins = [24]
    range_ = [[-180., 180.]]
    data = np.asarray(data)
    data, _ = np.histogramdd(data.T, range=range_, bins=bins)
    data = data / np.sum(data)
    return data


def get_distributions(data: dict) -> dict:

    # hist data is organized with resnames as keys and a tuple of histograms as values
    hist_data = {}
    # loop over residue names
    for i in data[0].keys():
        data1 = make_hist(data[0][i])
        data2 = make_hist(data[1][i])

        hist_data[i] = (data1, data2)

    return hist_data


def get_kbl(data: dict) -> dict:
    kbl_dict = {}

    for i in data:
        a = data[i][0]
        b = data[i][1]

        kbl = []

        for indx, item in np.ndenumerate(a):
            j = a[indx]
            k = b[indx]
            total = (j + k) / 2

            if j == 0.:
                j_term = 0.
            else:
                j_term = j * np.log(j / (total))

            if k == 0.:
                k_term = 0.
            else:
                k_term = k * np.log(k / (total))

            kbl.append((j_term + k_term) / 2)

        kbl = sum(kbl)
        kbl_dict[i] = kbl

    return kbl_dict
