#!/usr/bin/ env python3
# -*- coding: utf-8; fill-column: 100 truncate-lines: t -*-
#
# Trajectory simulation for molecular beam deflection by inhomogeneous electric fields
# Tool to plot spatial profiles of single quantum-state deflection simulations with CMIfly
#
# Copyright (C) 2015 Daniel Horke
# Copyright (C) 2020 Controlled Molecule Imaging Group, Center for Free-Electron Laser Science,
#                    Deutsches Elektronen-Synchrotron DESY and Universität Hamburg, Hamburg, Germany
#
# CMIfly is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, considering the amendment provided in LICENSE.md.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not,
# see <https://www.gnu.org/licenses/>.


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
parser.add_argument(
                    "-t",
                    "--traj",
                    action='store_true',
                    default='store_true',
                    help="plot trajectories")
args=parser.parse_args()

#define input path

PathIN = ''#/DataCloud/WorkCloud/Calculations/Dipeptide/'
#create correct systematic name for loadfile
qs = args.quantum_state.split(',')
loadfile = PathIN+args.molecule + "_" + str(args.voltage) + "kV.5fly"
nodepath = "_" + qs[0] + "_"+ qs[1] + "_"+ qs[2] + "_"+ qs[3] + "_" + str(args.isomer)

#open 5fly file and read data from correct node
f = tables.open_file(loadfile,'r')
simdata_array = f.get_node("/" + nodepath + "/FlightData")
xpos_final = [0 for x in range(len(simdata_array))]
ypos_final = [0 for x in range(len(simdata_array))]
time_final = [0 for x in range(len(simdata_array))]
xpos_inital = [0 for x in range(len(simdata_array))]
ypos_inital = [0 for x in range(len(simdata_array))]
time_inital = [0 for x in range(len(simdata_array))]

for i in range(len(simdata_array)):
    time_final[i] = simdata_array[i][0]
    xpos_final[i] = simdata_array[i][1]
    ypos_final[i] = simdata_array[i][2]
    time_inital[i] = simdata_array[i][8]
    xpos_inital[i] = simdata_array[i][9]
    ypos_inital[i] = simdata_array[i][10]


if args.profile !=None or args.save != None:
    hist, xedges = numpy.histogram(ypos_final, range=(-0.002,0.002),bins=args.bins)
    yposition = numpy.linspace(xedges[0]*1000, xedges[-1]*1000, num=args.bins)
    histx, xedgesx = numpy.histogram(xpos_final, range=(-0.002,0.002),bins=args.bins)
    xposition = numpy.linspace(xedgesx[0]*1000, xedges[-1]*1000, num=args.bins)

if args.scatter != None:
    plt.figure()
    plt.scatter(xpos_final,numpy.flipud(ypos_final))
    plt.show()

if args.profile != None:
    plt.figure()
    plt.plot(xposition,histx,marker='p')
    plt.xlim(xedgesx[0]*1000, xedgesx[-1]*1000)
    plt.xlabel('x-position (mm)')
    plt.ylabel('hits')
    plt.title(str(qs[0]) + "," + str(qs[1]) + "," + str(qs[2]) + "," + str(qs[3]) + "," + str(args.isomer) + " state of " + str(args.molecule) + " at " + str(args.voltage) + "kV")
    plt.show()

if args.traj != None:
    plt.figure()
    plt.plot(time_inital,xpos_inital,marker='*', label="inital postion")
    plt.plot(time_final,xpos_final,marker='p', label="final postion")
    plt.ylabel('x-position (mm)')
    plt.xlabel('time')
    plt.legend()
    plt.figure()
    plt.plot(xpos_inital,ypos_inital,marker='*', label="start postion")
    plt.plot(xpos_final,ypos_final,marker='p', label="stop postion")
    plt.ylabel('y-position (mm)')
    plt.xlabel('x-position (mm)')
    plt.legend()
    plt.show()

if args.save != None:
    savename = args.molecule + "_" + str(args.voltage) + "kV_" + qs[0] + "_"+ qs[1] + "_"+ qs[2] + "_"+ qs[3] + "_" + str(args.isomer) + ".dat"
    with open(savename, 'w') as f:
        for f1, f2 in zip(yposition, hist):
            print(f1,"\t",f2, file=f)
