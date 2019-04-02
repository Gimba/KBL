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
import numpy as np

from mpl_toolkits.mplot3d import axes3d, Axes3D
from hist_3D_lib import make_3dhist


def normalize_hist(hist, lines):
    hist = list(np.asarray(hist) / lines)
    return hist


def main(args):
    parser = argparse.ArgumentParser(description='Make a plot of how a distribution changes over time. Can be used to how distributions change over time.')
    parser.add_argument('file', help='input file in .xvg format')
    parser.add_argument('line_increment', help='increase the number of lines to read in by this parameter')
    parser.add_argument('increments', help='how often the number of lines should get bigger')
    parser.add_argument('bins', help='bins of the histogram')
    parser.add_argument('-final', dest='final', action='store_true', help='show histogram of all the lines')
    parser.add_argument('-hist3d', dest='hist3d', action='store_true', help='make 3d histogram, if not set make 2d '
                                                                            'histogram', default=False)
    parser.add_argument('-show', dest='show', action='store_true', help='show plot instead of saving to file',
                        default=False)
    args = parser.parse_args()

    line_increment = int(args.line_increment)
    increments = int(args.increments)
    n_bins = int(args.bins)

    with open(args.file, 'r') as f:

        data = []
        for line in f.readlines():
            # skip commented lines
            if '@' in line or '#' in line:
                continue
            data.append(float(line.split()[1]))

        all_hist_data = []

        title = args.file.strip('.xvg').split("/")[-1] + " angle distributions"

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
            n_frames = line_increment * increments + 1
            plt.hist2d(np.arange(0, n_frames), data[:n_frames], bins=[n_bins, 180])
            plt.title(title)
            plt.xlabel('frames')
            plt.ylabel('angle')

        if args.show:
            plt.show()
        else:
            plt.savefig(args.file.replace('.xvg', '.png'))



if __name__ == '__main__':
    main(sys.argv)
