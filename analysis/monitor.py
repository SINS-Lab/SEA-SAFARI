import os
import math
import time
import numpy as np
import platform
import argparse
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import safari_input
import subprocess
import detect_processor as detect
import tailer
import multiprocessing

#Some global variables
img = None
data = None
data_file = None
#This is needed as global otherwise it gets GC'd
ani = None
resX = 256
resY = 256
xi = 0
yi = 0
x_min = 0
x_max = 0
y_min = 0
y_max = 0

rx = 0
ry = 0

sx = 1
sy = 1

def updateplot(*args):
    q = args[1]
    global img
    global data
    try:
        changed = False
        while q.qsize() > 0:
            result = q.get_nowait()
            if result != '':
                x = result[0]
                y = result[1]
                data[y][x][0] = 1
                data[y][x][1] = 0
                data[y][x][2] = 0
                data[y][x][3] = 1
                changed = True
        if changed:
            img.set_array(data)
    except:
        pass

def plot(crystal, q):
    global img
    global data
    global ani
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    patches = []
    colours = []
    z_threshold = -safio.BDIST*1e3
    crystal = detect.sort_basis(crystal)
    for site in crystal:
        if site[2] < z_threshold:
            continue
        colours.append(site[2])
        circle = Circle((site[0]*sx, site[1]*sy), sx)
        patches.append(circle)

    p = PatchCollection(patches, cmap=plt.get_cmap('BuGn'), zorder = 5)
    p.set_array(np.array(colours))
    
    #Draw the basis
    ax.add_collection(p)
    #Add a heightmap
    fig.colorbar(p, ax=ax)

    norm = matplotlib.colors.Normalize(vmax=10, vmin=0, clip=True)
    img = plt.imshow(data, norm=norm, cmap=plt.get_cmap('plasma'), animated=True, zorder = 10)

    ani = anim.FuncAnimation(fig, updateplot, fargs=(q,), interval=50)
    plt.show()

def check_file(q):
    global data_file
    global data  
    for line in tailer.follow(data_file, delay=0.1):
        new_data = line.split()
        if len(new_data) > 3 and new_data[0] != 'X0':
            x = float(new_data[0])
            y = float(new_data[1])
            valid, x, y = to_data(x, y)
            if valid:
                q.put((x,y))

def to_data(x, y):
    x = resX * (x - x_min) / rx
    y = resY * (y - y_min) / ry
    if x >= 0 and x < resX and y >= 0 and y < resY:
        x = int(x)
        y = int(y)
        return True, x, y
    return False, 0, 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Data file to run from")
    parser.add_argument("-c", "--crys", help="Crystal file to use for lattice")
    parser.add_argument("-i", "--input", help="Input file to use for bounds")
    parser.add_argument("-s", "--size", help="Angular size of spot detector")
    parser.add_argument("-t", "--theta", help="Theta angle for detector")
    parser.add_argument("-e", "--emin", help="Minimum energy to consider")
    parser.add_argument("-r", "--emin_rel", help="Relative Minimum energy to consider")
    args = parser.parse_args()

    crysname = None
    inputname = None
    dataname = None

    #Check if we need to use -c and -i
    if not args.file.endswith('.data'):
        #Lets use this for crys, input and data
        crysname = args.file
        inputname = args.file + '.input'
        dataname = args.file + '.data'

    if crysname is not None and inputname is not None and dataname is not None:
        crystal = detect.loadCrystal(crysname)
        safio = safari_input.SafariInput(inputname)
        data = np.zeros((resY,resX, 4), dtype='float32')
        x_min = safio.XSTART
        x_max = safio.XSTOP
        y_min = safio.YSTART
        y_max = safio.YSTOP

        rx = x_max - x_min
        ry = y_max - y_min

        sx = resX / rx
        sy = resY / ry
        data_file = open(dataname)

        # Start the queue for handling multiprocessing
        q = multiprocessing.Queue()

        read_file = multiprocessing.Process(None, check_file, args=(q,))
        read_file.start()

        plot(crystal, q)
        input('')