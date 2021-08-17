import numpy as np  # Array manipulation
import argparse     # Argument parsing
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib   # Plotting
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt                    # Plotting
import matplotlib.image as mpimg                   # Plotting as image
import matplotlib.animation as anim                # Updating said image
from matplotlib.patches import Circle              # Circles on plot
from matplotlib.collections import PatchCollection # Also for said circles
import safari_input                                # Reading input files
import detect_processor as detect                  # Processing detections
import tailer                                      # Reading data file
import multiprocessing                             # Reading async from plotting

#Some global variables
img = None
impact_data = None  # Stores impact plot detections
e_theta_data = None # Stores E_Theta plot detections
data_file = None
#This is needed as global otherwise it gets GC'd
ani = None

#Resolution of the image
resX = 256
resY = 256

#Ranges for the plot
x_min = 0
x_max = 0
y_min = 0
y_max = 0
rx = 0
ry = 0

#Scaling factors for the plot
sx = 1
sy = 1

#Detector ranges
phi_min = -180
phi_max = 180
theta_min = 0
theta_max = 90
e_min = 0.5
e_max = 1e5

def update_impact_plot(*args):
    q = args[1]
    global img
    global impact_data
    try:
        changed = False
        while q.qsize() > 0:
            result = q.get_nowait()
            if result != '':
                x = result[0]
                y = result[1]
                #TODO make this lookup from colour map or something?

                if result[5]:
                    # If was marked yellow, reset it.
                    if impact_data[y][x][1] != 0:
                        impact_data[y][x][0] = 0
                        impact_data[y][x][1] = 0
                    # Add some redness
                    impact_data[y][x][0] += 0.25
                # Only mark it as yellow if not already red.
                elif impact_data[y][x][0] == 0:
                    impact_data[y][x][0] = 1
                    impact_data[y][x][1] = 1
                if impact_data[y][x][0] > 1:
                    impact_data[y][x][0] = 1
                if impact_data[y][x][1] > 1:
                    impact_data[y][x][1] = 1

                #The main reason we use this instead is for default transparency
                impact_data[y][x][3] = 0.75
                changed = True
        if changed:
            img.set_array(impact_data)
    except:
        pass

def impact_plot(crystal, q):
    global img
    global impact_data
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

    img = plt.imshow(impact_data, animated=True, zorder = 10)

    ani = anim.FuncAnimation(fig, update_impact_plot, fargs=(q,), interval=50)
    plt.show()

def check_file(q):
    global data_file
    global impact_data  
    for line in tailer.follow(data_file, delay=0.1):
        new_data = line.split()
        if len(new_data) > 5 and new_data[0] != 'X0':
            x = float(new_data[0])
            y = float(new_data[1])
            valid, x, y, E, theta, phi, detected = to_data(x, y, new_data)
            if valid:
                q.put((x,y, E, theta, phi, detected))

def to_data(x, y, data):
    x = resX * (x - x_min) / rx
    y = resY * (y - y_min) / ry
    #Check if it fits in the range
    if x >= 0 and x < resX and y >= 0 and y < resY:
        #Check if it fits in the detector
        E = float(data[3])
        theta = float(data[4])
        phi = float(data[5])
        detected =  E > e_min and E < e_max and \
                    theta > theta_min and theta < theta_max and \
                    phi > phi_min and phi < phi_max
        x = int(x)
        y = int(y)
        return True, x, y, E, theta, phi, detected
    return False, 0, 0, 0, 0, 0, False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Data file to run from")
    parser.add_argument("-c", "--crys", help="Crystal file to use for lattice")
    parser.add_argument("-i", "--input", help="Input file to use for bounds")
    parser.add_argument("-s", "--size", help="Angular size of spot detector")
    parser.add_argument("-t", "--theta", help="Theta angle for detector")
    parser.add_argument("-e", "--emin", help="Minimum energy to consider")
    args = parser.parse_args()

    crysname = None
    inputname = None
    dataname = None

    # Check if we need to use -c and -i
    if not args.file.endswith('.data'):
        # Lets use this for crys, input and data
        crysname = args.file
        inputname = args.file + '.input'
        dataname = args.file + '.data'

    if crysname is not None and inputname is not None and dataname is not None:
        crystal = detect.loadCrystal(crysname)
        safio = safari_input.SafariInput(inputname)
        impact_data = np.zeros((resY,resX, 4), dtype='float32')

        # Initialize the ranges
        x_min = safio.XSTART
        x_max = safio.XSTOP
        y_min = safio.YSTART
        y_max = safio.YSTOP

        rx = x_max - x_min
        ry = y_max - y_min

        sx = resX / rx
        sy = resY / ry

        asize = 10

        # Initialize angular ranges as well

        if args.size:
            asize = float(args.size)
        if args.theta:
            theta = float(args.theta)
            theta_min = theta - asize / 2.0
            theta_max = theta + asize / 2.0

        # This one is fixed by experimental geometry, so not configurable.
        phi_min = safio.PHI0 - asize / 2.0
        phi_max = safio.PHI0 + asize / 2.0


        data_file = open(dataname)

        # Start the queue for handling multiprocessing
        q = multiprocessing.Queue()

        read_file = multiprocessing.Process(None, check_file, args=(q,))
        read_file.start()

        impact_plot(crystal, q)
        input('')