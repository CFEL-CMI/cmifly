#!/opt/local/bin/python
# Tool to plot spatial profiles of single quantum-state deflection simulations with CMIfly
# Copyright (C) Daniel Horke 2015

import math, random, tarfile
import tables
from tables import *
import numpy
import argparse
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument(
                    "-v",
                    "--voltage",
                    type=float,
                    default=0.0,
                    help="Specify deflector voltage in kV. Default is 0.0kV")
parser.add_argument(
                    "-m",
                    "--molecule",
                    type=str,
                    default="OCS",
                    help="Specify molecule of interest. See code for implemented molecules.")
parser.add_argument(
                    "-q",
                    "--quantum_state",
                    default="0,0,0,0",
                    help="specify quantum state in format: J,Ka,Kc,M. Default is 0,0,0,0.")
parser.add_argument(
                    "-i",
                    "--isomer",
                    type=int,
                    default=0,
                    help="specify isomer to simulate. default = 0")
parser.add_argument(
                    "--scatter",
                    action='store_true',
                    default=None,
                    help="create scatter plot of final positions")
parser.add_argument(
                    "-p",
                    "--profile",
                    action='store_true',
                    default=None,
                    help="plot histogram in deflection direction")
parser.add_argument(
                    "-b",
                    "--bins",
                    type=int,
                    default=50,
                    help="number of bins for histogram. default = 50")
parser.add_argument(
                    "-s",
                    "--save",
                    action='store_true',
                    default=None,
                    help="save histogram data to ascii file")
args=parser.parse_args()

#define input path
cmiflypath = "" #'/afs/desy.de/user/h/horked/CMIfly/'

#create correct systematic name for loadfile
qs = args.quantum_state.split(',')
loadfile = args.molecule + "_" + str(args.voltage) + "kV.5fly"
nodepath = "_" + qs[0] + "_"+ qs[1] + "_"+ qs[2] + "_"+ qs[3] + "_" + str(args.isomer)

#open 5fly file and read data from correct node
f = tables.openFile(loadfile,'r')
simdata_array = f.getNode("/" + nodepath + "/FlightData")
xpos_final = [0 for x in range(len(simdata_array))]
ypos_final = [0 for x in range(len(simdata_array))]
for i in range(len(simdata_array)):
    xpos_final[i] = simdata_array[i][1]
    ypos_final[i] = simdata_array[i][2]

if args.profile !=None or args.save != None:
    hist, xedges = numpy.histogram(ypos_final, range=(-0.001,0.001),bins=args.bins)
    yposition = numpy.linspace(xedges[0]*1000, xedges[-1]*1000, num=args.bins)

if args.scatter != None:
    plt.figure(1)
    plt.scatter(xpos_final,ypos_final)
    plt.show()

if args.profile != None:
    plt.figure(2)
    plt.plot(yposition,hist,marker='p')
    plt.xlim(xedges[0]*1000, xedges[-1]*1000)
    plt.xlabel('y-position (mm)')
    plt.ylabel('hits')
    plt.title(str(qs[0]) + "," + str(qs[1]) + "," + str(qs[2]) + "," + str(qs[3]) + "," + str(args.isomer) + " state of " + str(args.molecule) + " at " + str(args.voltage) + "kV")
    plt.show()

if args.save != None:
    savename = args.molecule + "_" + str(args.voltage) + "kV_" + qs[0] + "_"+ qs[1] + "_"+ qs[2] + "_"+ qs[3] + "_" + str(args.isomer) + ".dat"
    with open(savename, 'w') as f:
        for f1, f2 in zip(yposition, hist):
            print(f1,"\t",f2, file=f)
