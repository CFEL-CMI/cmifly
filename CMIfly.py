#!/usr/bin/env python
# -*- coding: utf-8; truncate-lines: t -*-
# Trajectory simulation for molecular beam deflection by inhomogeneous electric fields
#
# This Python program is partly based on the peer-reviewed article
# Copyright (C) Yuan-Pin Chang, Daniel Horke, Sebastian Trippel and Jochen Küpper, CFEL, DESY 2014

# import generally necessary modules
import math, random, tarfile
import tables
import sys
from tables import *
import numpy as np
from scipy import interpolate
from scipy.integrate import ode
import argparse
import matplotlib.pyplot as plt

#-------------------------PARSE ARGUMENTS FROM COMMAND LINE
parser = argparse.ArgumentParser()
parser.add_argument(
                    "-v",
                    "--voltage",
                    type=float,
                    default=0.0,
                    help="Specify deflector voltage in kV. Default is 0.0kV")
parser.add_argument(
                    "-p",
                    "--particles",
                    type=int,
                    default=100,
                    help="Specify number of particles to simulate. Default is 100.")
parser.add_argument(
                    "-m",
                    "--molecule",
                    type=str,
                    default="adenine",
                    help="Specify molecule of interest. See code for implemented molecules.")
parser.add_argument(
                    "-s",
                    "--source",
                    type=str,
                    default="gauss",
                    help="Specify source samples. Gauss or uniform. default = gauss.")
parser.add_argument(
                    "-j",
                    "--jmax",
                    type=int,
                    default=0,
                    help="specify max J value to simulate. default = 0")
parser.add_argument(
                    "-i",
                    "--isomer",
                    type=int,
                    default=0,
                    help="specify isomer to simulate. default = 0")
parser.add_argument(
                    "-pf",
                    "--plotdeflectionfield",
                    action='store_true',
                    default=True,
                    help="plot deflection field")

args=parser.parse_args()
#-----------------------SETTING UP CALCULATION PARAMETERS
# number of molecules/trajectories to calculate
num_molecules = args.particles
# deflector voltage deifference between rod and trough (kV)
deflector_voltage = args.voltage
# Stark energy files
stark_filename = args.molecule + ".stark"
stark_path = "/usr/local/share/coldmol/molecules/"

outputfile = args.molecule + "_gauss_" + str(args.voltage) + "kV.5fly"
jmax=args.jmax

#-------------------------SETTING UP MACHINE PARAMETERS FOR NS-DYNAMIX
# molecular beam source position and spread [x, y, z, vx, yv, vz] (m, m, m, m/s, m/s, m/s)
# positions are relative to the beginning of the deflector at position (x, y, z) = (0, 0, 0)
# x- and y-width as well as vx- and vy-width must be identical
source_position = np.array([0, 0, -0.437, 0, 0, 670])
source_width = np.array([0.0001, 0.0001, 0, 1.5, 1.5, 7])
#source_width = np.array([0,0,0,0,0,0]) #for gauss sampling, these are st devs.

# deflector length (m)
deflector_length = 0.154
# electric field files, norm and gradient are saved separately
deflector_fieldnorm_filename     = 'deflector_field_norm'
deflector_fieldgradient_filename = 'deflector_field_gradient'
# deflector_voltage for which the fields were calculated (kV)
deflector_field_voltage = 5.
# geometric boundary of the deflector
rod_center = [0, 0.0036] # rod center (x,y)
rod_radius = 0.003 # rod radius
trough_center = [0, 0.002408] # trough center (x,y)
trough_radius = 0.0032 # trough radius
# detector_position (m)
detector_position = 0.301

# skimmer1
# skimmer1 position [x,y,z] and width
skimmer1_position = np.array([0,0,-0.362])
skimmer1_width = 0.002
# skimmer2
# skimmer2 position [x,y,z] and width
skimmer2_position = np.array([0,0,-0.02856])
skimmer2_width = 0.001
# skimmer3
# skimmer3 position [x,y,z] and width
skimmer3_position = np.array([0,0,0.167])
skimmer3_width = 0.0015

########## no input required below this line ##########

# convert units into internal units
# deflector voltage is only relevant as a scaling factor of the read fields
deflector_voltage_scaling = deflector_voltage / deflector_field_voltage
# mass (u -> kg)
u =  1.660538782e-27 # kg

#read molecule mass from stark file
f=tables.open_file(stark_path+stark_filename, 'r')
mass = f.get_node('/masses')[0][0]
f.close()
mass *= u

if args.source == "gauss":
    source_covariance = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i==j:
                source_covariance[i][j] = source_width[i]**2

def sample_source_uniform(position, width):
    """Sample the molecular beam source as defined above and return next molecules start position.

    This uses a random number generator to Monte-Carlo sample the defined source.
    """
    from numpy.random import uniform, normal
    pos = np.zeros(6)
    # we do geometric 2D uniform sampling in (x,y) to get a uniform 2d distribution
    while True:
        x, y = uniform(-0.5*width[0], 0.5*width[0], 2)
        if np.sqrt(x**2+y**2) <= 0.5*width[0]:
            pos[0:2] = x + position[0], y + position[1]
            break
    pos[2] = uniform(-0.5*width[2], 0.5*width[2]) + position[2]
    pos[3:6] = normal(position[3:6], 0.5*width[3:6])
    return pos

def sample_source_gauss(position, width):
    pos = np.zeros(6)
    pos = np.random.multivariate_normal(position,source_covariance)
    return pos

def derivative(t,position):
    """Calculate the derivative of the current phase-space position.

    d(x, y, z, vx, vy, vz) / dt = (vx, vy, vz, ax=dvx/dt, ay=dvy/dt, az=dvz/dt)

    For positions outside the deflector, ax = ay = az = 0.
    """
    acceleration_x, acceleration_y = acceleration(position) # obtain accerlation forces
    derivative = np.array(position)
    derivative[0:3] = position[3:6]
    #print('z:'+str(position[2]))
    #print('t:'+str(t))
    if position[2] >= 0.0 and position[2] <= deflector_length:
        ### this needs to be done right!
        derivative[3] = acceleration_x/mass
        derivative[4] = acceleration_y/mass
        derivative[5] = 0.0
    #print('acceleration:'+str(acceleration_x))
    else:
        derivative[3:6] = (0,0,0)
    #print('no acceleration')
    return derivative

def stark_effect(state, filename):
    """Read the starkeffect for |state| (J, Ka, Kc, M, isomer) from file |filename|

    Return
    Returns a 3-tuple of arrays containing the field strengths, the energies, and the effective dipole moments
    """
    f=tables.open_file(stark_path+filename, 'r')
    state_label =  "/_" + str(state[0]) + "/_" + str(state[1]) + "/_" + str(state[2]) + "/_" + str(state[3])  + "/_" + str(state[4])
    field_norm_array = f.get_node("/" + state_label + "/dcfield")
    energy_array = f.get_node("/" + state_label + "/dcstarkenergy")
    field_norm = np.array(field_norm_array.read())[0]
    energy = np.array(energy_array.read())[0]
    # calculate mueff from Stark energy
    mueff = np.zeros((len(field_norm),), np.float64)
    mueff[1:-1] = -1 * (energy[0:-2] - energy[2:]) / (field_norm[0:-2] - field_norm[2:])
    mueff[0] = 0.
    mueff[-1] = mueff[-2]
    mueff_interp = interpolate.interp1d(field_norm,mueff,bounds_error=False) # create interpolate object
    f.close()
    return mueff_interp

def read_deflection_field(filename):
    rawdata = open(filename,'r')
    # Read the file contents and generate a list with each line
    lines = rawdata.readlines()
    # initial positions of x,y,z grid
    line = lines[0][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xfs, yfs, zfs = map(float, linesplit)
    # step size for each dimension
    line = lines[1][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    dfx, dfy, dfz = map(float, linesplit)
    # number of steps for each dimension
    line = lines[2][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xfn, yfn, zfn = map(int, linesplit)
    # create a new array for import field grid values
    field = np.zeros((yfn,xfn),float)
    # create x,y mesh of field grid
    xf = np.linspace(xfs,xfs+dfx*(xfn-1),xfn)
    yf = np.linspace(yfs+dfy*(yfn-1),yfs,yfn) # transform -y to +y
    for xfi in range(xfn):
        for yfi in range(yfn):
            line = lines[3 + yfi + xfi * yfn][:-1] # remove EOL
            field[yfi,xfi] = float(line)
    rawdata.close()
    return xf, yf, field

def read_deflection_gradient(filename):
    f = open(filename,'r')
    # Read the file contents and generate a list with each line
    lines = f.readlines()
    # initial positions of x,y,z grid
    line = lines[0][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xgs, ygs, zgs = map(float, linesplit)
    # step size for each dimension
    line = lines[1][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    dxg, dyg, dzg = map(float, linesplit)
    # number of steps for each dimension
    line = lines[2][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xgn, ygn, zgn = map(int, linesplit)
    # create a new array for import field values
    gradx = np.zeros((ygn,xgn),float)
    grady = np.zeros((ygn,xgn),float)
    gradz = np.zeros((ygn,xgn),float)
    # create x,y mesh of field gridd
    xg = np.linspace(xgs,xgs+dxg*(xgn-1),xgn)
    yg = np.linspace(ygs+dyg*(ygn-1),ygs,ygn) # transform -y to +y
    for xgi in range(xgn):
        for ygi in range(ygn):
            line = lines[3 + ygi + xgi * ygn][:-2] # remove EOL and one space
            linesplit = line.split(" ")
            gradx[ygi][xgi], grady[ygi][xgi], gradz[ygi][xgi] = map(float, linesplit)
    return xg, yg, gradx, grady
    f.close()


def generate_acceleration_field(state, stark_filename, field_filename, gradient_filename, scaling=1.):
    """Generate a acceleration field at the same positions as the electric field is defined.

    Reads the Stark effect (energiers and effective dipole moments) from stark_filename.

    Reads the norm and gradient of the 2D electric field from field_filename and gradient_filename and (possibly) scales
    the field according to the scaling factor.

    Calculates a 2D acceleration field at the same positions as where the field is defined using equiation XXX from the
    manuscript.

    Returns generated interpolation objects of acceleration fields

    """
    # obtain field and gradient arrays
    field_x_grid, field_y_grid, deflection_field = read_deflection_field(field_filename)
    gradient_x_grid, gradient_y_grid, deflection_gradient_x, deflection_gradient_y = read_deflection_gradient(gradient_filename)
    # create interpolate objects
    deflection_field_interp = interpolate.interp2d(field_x_grid,field_y_grid,deflection_field*scaling,kind='linear')
    deflection_gradient_x_interp = interpolate.interp2d(gradient_x_grid,gradient_y_grid,deflection_gradient_x*scaling,kind='linear')
    deflection_gradient_y_interp = interpolate.interp2d(gradient_x_grid,gradient_y_grid,deflection_gradient_y*scaling,kind='linear')
    # and so forth
    return deflection_field_interp, deflection_gradient_x_interp, deflection_gradient_y_interp

def acceleration(position):
    """Calculate acceleration at a certain position in the current acceleration field."""
    deflection_gradient_x_value = deflection_gradient_x(position[0],position[1])[0]
    deflection_gradient_y_value = deflection_gradient_y(position[0],position[1])[0]
    mueff_value = mueff(deflection_field(position[0],position[1])[0])
    return deflection_gradient_x_value * mueff_value, deflection_gradient_y_value * mueff_value


def fly(initial_position, acceleration_z_bounds, skimmer1, skimmer2, skimmer3, final_z, final_t, norm_field_interp,norm_field_interp_0 ):
    """Propagate a molecule from it's initial phase-space |position| through the |acceleration| field (defined within its
    z_bounds) until the molecule reaches the final_z position ot final_t time is elapsed.
    """
    deflector_start, deflector_length = acceleration_z_bounds
    skimmer1_position, skimmer1_width = skimmer1
    skimmer2_position, skimmer2_width = skimmer2
    skimmer3_position, skimmer3_width = skimmer3
    dt = 1e-6
    t = 0. # assume time initial time is zero
    hit_deflector = False
    hit_skimmer1 = False
    hit_skimmer2 = False
    hit_skimmer3 = False
    flag=0
    integral=ode(derivative)
    integral.set_integrator('vode',method='BDF',with_jacobian=False,atol=1e-8,rtol=1e-8,nsteps=100000)
    integral.set_initial_value(initial_position,t)
    field_x_grid, field_y_grid, deflection_field_norm = read_deflection_field(deflector_fieldnorm_filename)
    field_limit = np.max(deflection_field_norm)/10000
    while integral.successful() and integral.y[2] <= final_z and integral.t<= final_t and not hit_deflector and not hit_skimmer1 and not hit_skimmer2 and not hit_skimmer3:
        integral.integrate(integral.t + dt)
        #print('calculated time:'+str(integral.t))
        if rod_radius != 0:
            hit_deflector = integral.y[2] > deflector_start and integral.y[2] < (deflector_length+deflector_start) and \
                        (rod_radius**2 > (rod_center[0]-integral.y[0])**2 + (rod_center[1]-integral.y[1])**2 or \
                        trough_radius**2 < (trough_center[0]-integral.y[0])**2 + (trough_center[1]-integral.y[1])**2)
        if rod_radius == 0 and deflector_voltage ==0:
            hit_deflector = integral.y[2] > deflector_start and integral.y[2] < (deflector_length+deflector_start) and \
                             norm_field_interp_0(integral.y[1],integral.y[0]) <= field_limit
        else:
            hit_deflector = integral.y[2] > deflector_start and integral.y[2] < (deflector_length+deflector_start) and \
                             norm_field_interp(integral.y[1],integral.y[0]) <= field_limit
        hit_skimmer1 = integral.y[2] > skimmer1_position[2] and integral.y[2] < (skimmer1_position[2]+0.003) and \
                      skimmer1_width**2/4 < (skimmer1_position[0]-integral.y[0])**2 + (skimmer1_position[1]-integral.y[1])**2
        hit_skimmer2 = integral.y[2] > skimmer2_position[2] and integral.y[2] < (skimmer2_position[2]+0.003) and \
                      skimmer2_width**2/4 < (skimmer2_position[0]-integral.y[0])**2 + (skimmer2_position[1]-integral.y[1])**2
        hit_skimmer3 = integral.y[2] > skimmer3_position[2] and integral.y[2] < (skimmer3_position[2]+0.003) and \
                       skimmer3_width**2/4 < (skimmer3_position[0]-integral.y[0])**2 + (skimmer3_position[1]-integral.y[1])**2
        if integral.y[2]>=0 and flag==0:
            integral.set_initial_value(integral.y,integral.t)
            flag=1
    #print(integral.y, integral.t)
    return integral.y, integral.t

#defining column headers for hdf5 file
class Particle(IsDescription):
    time_final        = FloatCol(pos=1)
    xpos_final        = FloatCol(pos=2)
    ypos_final        = FloatCol(pos=3)
    zpos_final        = FloatCol(pos=4)
    xvel_final        = FloatCol(pos=5)
    yvel_final        = FloatCol(pos=6)
    zvel_final        = FloatCol(pos=7)
    time_initial      = FloatCol(pos=8)
    xpos_initial      = FloatCol(pos=9)
    ypos_initial      = FloatCol(pos=10)
    zpos_initial      = FloatCol(pos=11)
    xvel_initial      = FloatCol(pos=12)
    yvel_initial      = FloatCol(pos=13)
    zvel_initial      = FloatCol(pos=14)

def plot_norm(field, x, y,lab):
    plt.figure()
    plt.pcolormesh(x*mm,y*mm,field/Vm_kVcm, cmap='plasma')
    plt.colorbar()
    plt.title(lab+'Field (kV/cm)')
    plt.show()

#print message to user
print("Runnig CMIfly for ", args.molecule)
print("Using a ",args.source," source distribution to fly ",args.particles," particles per state")
print("Calculating up to J state ",args.jmax, '\n')

if "__main__" == __name__:
    output = open_file(outputfile, mode = "w", title = (str(args.molecule) + " at " + str(args.voltage) + "kV"))
    for countj in range(0,jmax+1):
        for countka in range(0,countj+1):
            for countkc in range(0,countj+1):
                for countm in range(0,countj+1):
                    if countj !=0 and countka+countkc == 0:
                        pass
                    elif countka+countkc > countj+1 or countka+countkc < countj:
                        pass
                    else:
                        print("Currently flying for state ",countj, countka, countkc, countm)
                        group = output.create_group("/", ("_"+str(countj)+"_"+str(countka)+"_"+str(countkc)+"_"+str(countm)+"_"+str(args.isomer)), ("J="+str(countj)+" Ka="+str(countka)+" Kc="+str(countkc)+" M="+str(countm)+" Isomer="+str(args.isomer)))
                        quantum_state =(countj, countka, countkc, countm, args.isomer)
                        # generate acceleration file from Stark-effect file and quantum state, and deflection field norm and gradient
                        acceleration_fields = generate_acceleration_field(quantum_state, stark_filename,
                                                                   deflector_fieldnorm_filename, deflector_fieldgradient_filename,
                                                                   deflector_voltage_scaling)
                        deflection_field, deflection_gradient_x, deflection_gradient_y = acceleration_fields
                        acceleration_fields_0 = generate_acceleration_field(quantum_state, stark_filename,deflector_fieldnorm_filename, deflector_fieldgradient_filename,1)
                        deflection_field_0, deflection_gradient_x_0, deflection_gradient_y_0 = acceleration_fields_0
                        mueff = stark_effect(quantum_state, stark_filename)
                        # hit count
                        hit = 0
                        # open output file
                        table = output.create_table("/"+("_"+str(countj)+"_"+str(countka)+"_"+str(countkc)+"_"+str(countm)+"_"+str(args.isomer)), 'FlightData', Particle)
                        # let them fly
                        for n in range(num_molecules):
                            if args.source == "uniform":
                                initial = sample_source_uniform(source_position, source_width)
                            else:
                                initial = sample_source_gauss(source_position, source_width)
                            final, time = fly(initial, (0, deflector_length), (skimmer1_position, skimmer1_width), (skimmer2_position, skimmer2_width), (skimmer3_position, skimmer3_width),detector_position, 0.1, deflection_field, deflection_field_0)
                            if final[2] >= detector_position:
                                # trajectory reached end of deflector -> write initial and final positions to output file
                                # save final phasespace in output file
                                table.append([(str(time), final[0],final[1],final[2],final[3],final[4],final[5],str(0.),initial[0],initial[1],initial[2],initial[3],initial[4],initial[5])])
                                table.flush()
                            else:
                                hit += 1
                            perc_progress = n/num_molecules * 100
                            if perc_progress.is_integer():
                                sys.stdout.write("Progress: " + str(np.round(perc_progress)) + "% \r")
                                sys.stdout.flush()
                        print("DONE, ", str(hit)," particles lost and ",str(num_molecules-hit)," particles detected! \n \n")
                        # all molecules flown, close output and finish
    output.close()

if args.plotdeflectionfield == True:
    field_x_grid, field_y_grid, deflection_field_norm = read_deflection_field(deflector_fieldnorm_filename)
    plot_norm(deflection_field(field_x_grid,field_y_grid), field_x_grid, field_y_grid, 'interpolated')
    field_x_grid_original, field_y_grid_original, deflection_field_original = read_deflection_field(deflector_fieldnorm_filename)
    plot_norm(deflection_field_original*deflector_voltage_scaling,field_x_grid_original, field_y_grid_original, 'original')
    step_x, step_y, field_x_grad, field_y_grad = read_deflection_gradient(deflector_fieldgradient_filename)
    plot_norm(field_x_grad,step_x, step_y, 'original')
    plot_norm(field_y_grad,step_x, step_y, 'original')
    plot_norm(deflection_gradient_x(field_x_grid,field_y_grid),step_x, step_y, 'interpolated')
    plot_norm(deflection_gradient_y(field_x_grid,field_y_grid),step_x, step_y, 'interpolated')
