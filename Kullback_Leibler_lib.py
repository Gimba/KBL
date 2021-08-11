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

import numpy as np
import pytraj as pt
from scipy.spatial.distance import jensenshannon


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


def extract_angles(files, residues, angles, topology):
    """

    :param files: trajectories
    :param residues: range of residues given in cpptraj convention
    :param angles: psi,phi,chi by default
    :param topology
    :return: data_dict: contains angles with residue identifier as keys and values as list
    """
    # "-1" needed for conversion between python convention (starts with 0) and cpptraj convention (starts with 1)
    resids = [int(r[3:]) - 1 for r in residues]

    # construct this to map resids back to full res name
    residue_dict = {}
    for r in residues:
        residue_dict[int(r[3:])] = r

    traj = pt.TrajectoryIterator(files, topology)
    data = pt.multidihedral(traj, dihedral_types=' '.join(angles), resrange=resids)
    # convert data to dictionary
    data_dict = {}
    for k in data.keys():
        data_dict[residue_dict[int(k.split(":")[1])] + " " + k.split(":")[0].replace("chip", "chi")] = np.asarray(
            data[k].values)

    return data_dict


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


def read_in_data(files1: "", files2: "", mutations: list = [], angles: list = [], residues: list = [],
                 topologies: list = []) -> dict:
    angle_mutual_residues = get_resids_from_files(files1, files1, angles, residues)

    data = {}

    def add_mutation_res():
        angles_1 = np.array(read_1_set(files1, end1, mutations[0], angles)).T
        angles_2 = np.array(read_1_set(files2, end2, mutations[1], angles)).T

        # handle case where the mutation has not got the chi1 angle
        if len(angles_1) > len(angles_2):
            angles_1 = np.delete(angles_1, -1, 0)
        # handle case where the non-mutation has not got the chi1 angle
        elif len(angles_1) < len(angles_2):
            angles_2 = np.delete(angles_2, -1, 0)
        # naming according to reference residue ids (dir1 pdb)
        data[mutations[0]] = (angles_1, angles_2)

    # # add specified mutation residue which has different file names
    if len(mutations) > 1:
        if residues:
            # check if given mutation is in specified residues
            if int(mutations[0][3:]) in residues:
                add_mutation_res()

        else:
            add_mutation_res()

    # for res in angle_mutual_residues[angles[0]]:
    dir1_angles = extract_angles(files1, angle_mutual_residues[angles[0]], angles, topologies[0])

    dir2_angles = extract_angles(files2, angle_mutual_residues[angles[0]], angles, topologies[1])
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


def get_jsd(data: dict) -> dict:
    jsd_dict = {}

    for i in data:
        a = data[i][0]
        b = data[i][1]
        jsd = jensenshannon(a, b)
        # kbl = []
        #
        # for indx, item in np.ndenumerate(a):
        #     j = a[indx]
        #     k = b[indx]
        #     total = (j + k) / 2
        #
        #     if j == 0.:
        #         j_term = 0.
        #     else:
        #         j_term = j * np.log(j / (total))
        #
        #     if k == 0.:
        #         k_term = 0.
        #     else:
        #         k_term = k * np.log(k / (total))
        #
        #     kbl.append((j_term + k_term) / 2)
        #
        # kbl = sum(kbl)
        # print(kbl)
        jsd_dict[i] = jsd
    return jsd_dict
