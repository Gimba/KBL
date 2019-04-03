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
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from mpl_toolkits.mplot3d import axes3d, Axes3D
from hist_3D_lib import make_3dhist

color_map = 'viridis'
plt.rcParams['axes.facecolor'] = cm.get_cmap(color_map)(0)

def normalize_hist(hist, lines):
    hist = list(np.asarray(hist) / lines)
    return hist


def main(args):
    parser = argparse.ArgumentParser(description='Make a plot of how a distribution changes over time. Can be used to how distributions change over time.')
    parser.add_argument('-f', dest='file', help='input file in .xvg format')
    # TODO: change bins to ybins and li,i to xbins
    parser.add_argument('-b', dest='bins', help='bins of the histogram', default=False)
    parser.add_argument('-li', dest='line_increment', nargs='?', help='increase the number of lines to read in by '
                                                                      'this parameter', default=False)
    parser.add_argument('-i', dest='increments', nargs='?', help='how often the number of lines should get bigger',
                        default=False)
    parser.add_argument('-final', dest='final', action='store_true', help='show histogram of all the lines')
    parser.add_argument('-hist3d', dest='hist3d', action='store_true', help='make 3d histogram, if not set make 2d '
                                                                            'histogram', default=False)
    parser.add_argument('-show', dest='show', action='store_true', help='show plot instead of saving to file',
                        default=False)
    parser.add_argument('-pre', nargs='?', help='text that appears before the title', default="")
    args = parser.parse_args()



    out_file = args.file.replace('.xvg', '.png')

    title = args.pre + args.file.strip('.xvg').split("/")[-1]

    with open(args.file, 'r') as f:

        data = []
        for line in f.readlines():
            # skip commented lines
            if '@' in line or '#' in line:
                continue
            data.append(float(line.split()[1]))

        if not (args.bins or args.increments or args.line_increment):
            increments = 10
            line_increment = len(data) // increments
            n_bins = 30
            print('using default settings for parameters bins (30), increments (10) and calculated line_increment (' +
                  str(line_increment) + ')')

        else:
            line_increment = int(args.line_increment)
            increments = int(args.increments)
            n_bins = int(args.bins)

        if len(data) < increments * line_increment:
            line_increment = len(data) // increments
        all_hist_data = []

        if args.hist3d:
            bins = tuple([n_bins])
            range_ = tuple([tuple([-180., 180.])])
            frames_counts = []
            for i in range(1, increments + 1):
                n_frames = line_increment * i
                frames_counts.append(n_frames)
                hist, _ = np.histogramdd([data[:n_frames]], range=range_, bins=bins)
                hist = normalize_hist(hist, line_increment * i)
                all_hist_data += list(hist)
            if args.final:
                hist, _ = np.histogramdd([data], range=range_, bins=bins)
                hist = normalize_hist(hist, len(data))
                all_hist_data += list(hist)
                increments += 1
                frames_counts.append("final")

            make_3dhist(n_bins, increments, all_hist_data, frames_counts, title)


        else:
            n_frames = line_increment * increments
            plt.hist2d(np.arange(0, n_frames), data[:n_frames], bins=[n_bins, 180], cmap=color_map)
            plt.ylim(-180, 180)
            plt.title(title)
            plt.xlabel('frames')
            plt.ylabel('angle')

        if args.show:
            plt.show()
        else:
            if args.hist3d:
                out_file = out_file.split('/')
                out_file = '/'.join(out_file[:-1]) + '/hist3d.' + out_file[-1]
            plt.savefig(out_file)
            print('saved plot to:', out_file)



if __name__ == '__main__':
    main(sys.argv)
