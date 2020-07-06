#!/usr/bin/ env python3
# -*- coding: utf-8; fill-column: 100 truncate-lines: t -*-
#
# Trajectory simulation for molecular beam deflection by inhomogeneous electric fields
# Programm for evaluating thermal population of J states.
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


import tables
import numpy
import math
from tables import *
import argparse
import matplotlib
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument(
                    "-t",
                    "--temperature",
                    type=float,
                    default=1.0,
                    help="Specify temperature of thermal population analysis. Default = 1.0 K.")
parser.add_argument(
                    "-m",
                    "--molecule",
                    type=str,
                    default="OCS",
                    help="Specify molecule of interest. See code for implemented molecules.")
parser.add_argument(
                    "-v",
                    "--voltage",
                    type=float,
                    default=0.0,
                    help="Specify deflector voltage in kV. Default is 0.0kV")
parser.add_argument(
                    "-j",
                    "--jmax",
                    type=int,
                    default=10,
                    help="specify max J value. default = 10")
parser.add_argument(
                    "-i",
                    "--isomer",
                    type=int,
                    default=0,
                    help="specify isomer to simulate. default = 0")
parser.add_argument(
                    "--population",
                    action = 'store_true',
                    default = None,
                    help="print population analysis on screen if called.")
parser.add_argument(
                    "-p",
                    "--plot",
                    action = 'store_true',
                    default = None,
                    help="plot population analysis if called.")
parser.add_argument(
                    "-b",
                    "--bins",
                    type=int,
                    default=50,
                    help="number of bins for histogram. default = 50")
parser.add_argument(
                    "-r",
                    "--range",
                    default="-1,1",
                    help="enter range of histogram in mm in format x,y. default is -1,1.")
parser.add_argument(
                    "-s",
                    "--save",
                    action = 'store_true',
                    default = None,
                    help="save population analysis if called.")
args=parser.parse_args()

boltzmann = 1.38065053e-23

#position of the stark file
stark_filename = args.molecule + ".stark"
stark_dir = ""#/usr/local/share/coldmol/molecules/"
flyfile = args.molecule + "_" +"gauss_" + str(args.voltage) + "kV.5fly"
fstark=tables.open_file(stark_dir+stark_filename, 'r') # stark files
ffly = tables.open_file(flyfile, 'r')

#histogram parameters
hist_range = args.range.split(",")
for i in range(2): hist_range[i] = int(hist_range[i])/1000
profile_scaled = [0 for x in range(args.bins)]

total_population = [0 for x in range(args.jmax+1)]
fractional_population = [0 for x in range(args.jmax+1)]
population_sum = 0
for J in range(0, args.jmax+1):
    population_J = 0;
    for Ka in range(0,J+1):
        for Kc in range(0,J+1):
            for M in range(0,J+1):
                if J !=0 and Ka+Kc == 0:
                    pass
                elif Ka+Kc > J+1 or Ka+Kc < J:
                    pass
                else:
                    array=fstark.get_node("/_"+str(J)+"/_"+str(Ka)+"/_"+str(Kc)+"/_"+str(M)+"/_0/dcstarkenergy")
                    starkcurve=numpy.array(array.read())
                    simdata_array = ffly.get_node("/" + ( "_" + str(J) + "_"+ str(Ka) + "_"+ str(Kc) + "_"+ str(M) + "_" + str(args.isomer)) + "/FlightData")
                    ypos_final = [0 for x in range(len(simdata_array))]
                    for i in range(len(simdata_array)):
                        # the deflection is in the x-direction, so ypos_final should be asigned to
                        # the column that contains the xpos_final
                        ypos_final[i] = simdata_array[i][1]
                    profile, xedges = numpy.histogram(ypos_final, range=hist_range ,bins=args.bins)
                    energy = starkcurve[0][0]
                    # Nuclear spin statistical weights should be checked for each molecule
                    nssw = 1.0
                    if M == 0:
                        Mdegen = 1.0
                    else:
                        Mdegen = 2.0
                    population = nssw * Mdegen * math.exp(-energy/(boltzmann*args.temperature));
                    for i in range(args.bins):
                        profile_scaled[i] = profile_scaled[i] + (profile[i] * population)
                    population_J = population_J + population
    population_sum = population_sum + population_J
    total_population[J] = population_J
for i in range(args.jmax+1):
   fracpop = total_population[i]/population_sum
   fractional_population[i] = fracpop

if args.plot != None:
    plt.figure(1)
    yposition = numpy.linspace(xedges[0]*1000, xedges[-1]*1000, num=args.bins)
    plt.plot(yposition, profile_scaled,marker='o')
    plt.show()

if args.save != None:
    yposition = numpy.linspace(xedges[0]*1000, xedges[-1]*1000, num=args.bins)
    savename = args.molecule + "_" + str(args.voltage) + "kV_" + str(args.temperature) + "K.dat"
    savefile = open(savename, 'w')
    for i in range(args.bins):
        savefile.write(str(yposition[i]) + " " + str(profile_scaled[i]) + "\n")
    savefile.close()

if args.population != None:
    template = "{0:10}{1:25}{2:25}"
    print(("Thermal population analysis of " + str(args.molecule) + " at " + str(args.temperature) + "K"))
    print((template.format("J", "Total Population", "Fractional Population")))
    for i in range(args.jmax+1):
        print(template.format(str(i), str(total_population[i]),str(fractional_population[i])))
