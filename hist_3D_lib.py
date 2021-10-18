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


import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['axes.formatter.useoffset'] = False
def make_3dhist(n_bins, n_histograms, data, frame_counts, title):
    fig = plt.figure(figsize=(plt.figaspect(0.5))*2)



    #===============
    #  First subplot
    #===============
    # set up the axes for the first plot
    ax = fig.add_subplot(1, 2, 1, projection='3d')

    # intervals = 360//n_bins
    # loc = plticker.MultipleLocator(base=intervals)
    # # ax.xaxis.set_major_locator(loc)
    # ax.yaxis.set_major_locator(loc)
    # ax.grid(which='major', axis='both', linestyle='-')

    dx = (1 / 3) * n_bins
    dy = (1 / 36) * n_histograms
    dz = data  # my data

    xpos = (list(range(-180,180, 360 //n_bins)))*n_histograms
    xpos = list(np.asarray(xpos) + (dx))
    temp = []

    ypos = []
    for i in range(1, n_histograms+1):
        temp = [i]*n_bins
        ypos+=temp
    ypos = list(np.asarray(ypos) - dy*1.5)
    zpos = np.zeros_like(xpos)

    cmap = cm.get_cmap('rainbow')  # Get desired colormap - you can change this!
    max_height = 360  # get range of colorbars so we can normalize
    min_height = 0
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k*n_bins - min_height) / max_height) for k in range(n_bins)]

    rgba = rgba*n_histograms



    ax.bar3d(xpos, ypos, zpos, dx, dy, dz,color=rgba)
    plt.yticks(np.arange(1,n_histograms+1,1),frame_counts, verticalalignment='baseline',
                   horizontalalignment='left', rotation=-15,)
    plt.xticks(range(-180,181,n_bins), range(-180,181,n_bins))
    plt.title(title)
    plt.xlabel("angle")
    plt.ylabel("frames",labelpad=10)
    ax.set_facecolor('white')

    ax.set_zlabel("%")
    # plt.show()
    if type(frame_counts[-1]) == int:
        max_frames = str(frame_counts[-1])
    else:
        max_frames = str(frame_counts[-2])

    plt.savefig(title.split()[0] + "_" + max_frames + ".png",  bbox_inches='tight', pad_inches=0)
    return