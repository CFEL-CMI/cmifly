#!/opt/local/bin/python
# Trajectory simulation for molecular beam deflection by inhomogeneous electric fields
#
# This Python program is partly based on the peer-reviewed article
# Copyright (C) Yuan-Pin Chang, Daniel Horke, Sebastian Trippel and Jochen Küpper, CFEL, DESY 2014
# all further developments (C) Daniel Horke, 2015

# import generally necessary modules
import math, random, tarfile
import tables
from tables import *
import numpy
from scipy import interpolate
from scipy.integrate import odeint
import argparse

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
                    default="OCS",
                    help="Specify molecule of interest. See code for implemented molecules.")
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
args=parser.parse_args()
#-----------------------SETTING UP CALCULATION PARAMETERS
# number of molecules/trajectories to calculate
num_molecules = args.particles
# deflector voltage deifference between rod and trough (kV)
deflector_voltage = args.voltage
# Stark energy files
stark_filename = args.molecule + ".stark"
# mass of molecule (u -- unified atomic mass units)
if args.molecule == "OCS":
    mass = 60.075
elif args.molecule == "mephenesin":
    mass = 182.216
elif args.molecule == "water":
    mass = 18.015
else:
    print("unrecognised molecule!")
outputfile = args.molecule + "_" + str(args.voltage) + "kV.5fly"
jmax=args.jmax

#-------------------------SETTING UP MACHINE PARAMETERS FOR NS-DYNAMIX
# molecular beam source position and spread [x, y, z, vx, yv, vz] (m, m, m, m/s, m/s, m/s)
# positions are relative to the beginning of the deflector at position (x, y, z) = (0, 0, 0)
# x- and y-width as well as vx- and vy-width must be identical
source_position = numpy.array([0, 0, -0.2462, 0, 0, 1900])
source_width = numpy.array([0.001, 0.001, 0, 1, 1, 19])

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
detector_position = 0.303

# skimmer1
# skimmer1 position [x,y,z] and width
skimmer1_position = numpy.array([0,0,-0.193])
skimmer1_width = 0.002
# skimmer2
# skimmer2 position [x,y,z] and width
skimmer2_position = numpy.array([0,0,-0.02856])
skimmer2_width = 0.001
# skimmer3
# skimmer3 position [x,y,z] and width
skimmer3_position = numpy.array([0,0,0.167])
skimmer3_width = 0.0015

########## no input required below this line ##########

# test input for consistency
assert(source_width[0] == source_width[1])
assert(source_width[3] == source_width[4])

# convert units into internal units
# deflector voltage is only relevant as a scaling factor of the read fields
deflector_voltage_scaling = deflector_voltage / deflector_field_voltage
# mass (u -> kg)
u =  1.660538782e-27 # kg
mass *= u

def sample_source(position, width):
    """Sample the molecular beam source as defined above and return next molecules start position.

    This uses a random number generator to Monte-Carlo sample the defined source.
    """
    from numpy.random import uniform, normal
    pos = numpy.zeros(6)
    # we do geometric 2D uniform sampling in (x,y) to get a uniform 2d distribution
    while True:
        x, y = uniform(-0.5*width[0], 0.5*width[0], 2)
        if numpy.sqrt(x**2+y**2) <= 0.5*width[0]:
            pos[0:2] = x + position[0], y + position[1]
            break
    pos[2] = uniform(-0.5*width[2], 0.5*width[2]) + position[2]
    pos[3:6] = normal(position[3:6], 0.5*width[3:6])
    return pos

def derivative(position, t):
    """Calculate the derivative of the current phase-space position.

    d(x, y, z, vx, vy, vz) / dt = (vx, vy, vz, ax=dvx/dt, ay=dvy/dt, az=dvz/dt)

    For positions outside the deflector, ax = ay = az = 0.
    """
    acceleration_x, acceleration_y = acceleration(position) # obtain accerlation forces
    derivative = numpy.array(position)
    derivative[0:3] = position[3:6]
    if position[2] >= 0.0 and position[2] <= deflector_length:
        ### this needs to be done right!
        derivative[3] = acceleration_x/mass
        derivative[4] = acceleration_y/mass
        derivative[5] = 0.0
    else:
        derivative[3:6] = (0,0,0)
    return derivative

def stark_effect(state, filename):
    """Read the starkeffect for |state| (J, Ka, Kc, M, isomer) from file |filename|

    Return
    Returns a 3-tuple of arrays containing the field strengths, the energies, and the effective dipole moments
    """
    f=tables.openFile(filename, 'r')
    state_label =  "/_" + str(state[0]) + "/_" + str(state[1]) + "/_" + str(state[2]) + "/_" + str(state[3])  + "/_" + str(state[4])
    field_norm_array = f.getNode("/" + state_label + "/dcfield")
    energy_array = f.getNode("/" + state_label + "/dcstarkenergy")
    field_norm = numpy.array(field_norm_array.read())[0]
    energy = numpy.array(energy_array.read())[0]
    # calculate mueff from Stark energy
    mueff = numpy.zeros((len(field_norm),), numpy.float64)
    mueff[1:-1] = -1 * (energy[0:-2] - energy[2:]) / (field_norm[0:-2] - field_norm[2:])
    mueff[0] = 0.
    mueff[-1] = mueff[-2]
    mueff_interp = interpolate.interp1d(field_norm,mueff) # create interpolate object
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
    field = numpy.zeros((xfn,yfn),float)
    # create x,y mesh of field grid
    xf = numpy.linspace(xfs,xfs+dfx*(xfn-1),xfn)
    yf = numpy.linspace(yfs+dfy*(yfn-1),yfs,yfn) # transform -y to +y
    for xfi in range(xfn):
        for yfi in range(yfn):
            line = lines[3 + yfi + xfi * yfn][:-1] # remove EOL
            field[xfi][yfi] = float(line)
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
    gradx = numpy.zeros((xgn,ygn),float)
    grady = numpy.zeros((xgn,ygn),float)
    gradz = numpy.zeros((xgn,ygn),float)
    # create x,y mesh of field gridd
    xg = numpy.linspace(xgs,xgs+dxg*(xgn-1),xgn)
    yg = numpy.linspace(ygs+dyg*(ygn-1),ygs,ygn) # transform -y to +y
    for xgi in range(xgn):
        for ygi in range(ygn):
            line = lines[3 + ygi + xgi * ygn][:-2] # remove EOL and one space
            linesplit = line.split(" ")
            gradx[xgi][ygi], grady[xgi][ygi], gradz[xgi][ygi] = map(float, linesplit)
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

def fly(initial_position, acceleration_z_bounds, skimmer1, skimmer2, skimmer3, final_z, final_t):
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
    position = odeint(derivative, initial_position, [t, t+dt])
    while position[-1,2] <= final_z and t + dt <= final_t and not hit_deflector and not hit_skimmer1 and not hit_skimmer2 and not hit_skimmer3:
        t += dt
        position = odeint(derivative, position[-1,:], [t, t+dt])
        hit_deflector = position[-1,2] > deflector_start and position[-1,2] < (deflector_length+deflector_start) and \
                        (rod_radius**2 > (rod_center[0]-position[-1,0])**2 + (rod_center[1]-position[-1,1])**2 or \
                        trough_radius**2 < (trough_center[0]-position[-1,0])**2 + (trough_center[1]-position[-1,1])**2)
        hit_skimmer1 = position[-1,2] > skimmer1_position[2] and position[-1,2] < (skimmer1_position[2]+0.003) and \
                      skimmer1_width**2/4 < (skimmer1_position[0]-position[-1,0])**2 + (skimmer1_position[1]-position[-1,1])**2
        hit_skimmer2 = position[-1,2] > skimmer2_position[2] and position[-1,2] < (skimmer2_position[2]+0.003) and \
                      skimmer2_width**2/4 < (skimmer2_position[0]-position[-1,0])**2 + (skimmer2_position[1]-position[-1,1])**2
        hit_skimmer3 = position[-1,2] > skimmer3_position[2] and position[-1,2] < (skimmer3_position[2]+0.003) and \
                       skimmer3_width**2/4 < (skimmer3_position[0]-position[-1,0])**2 + (skimmer3_position[1]-position[-1,1])**2
    return position[-1,:], t

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
                        group = output.create_group("/", ("_"+str(countj)+"_"+str(countka)+"_"+str(countkc)+"_"+str(countm)+"_"+str(args.isomer)), ("J="+str(countj)+" Ka="+str(countka)+" Kc="+str(countkc)+" M="+str(countm)+" Isomer="+str(args.isomer)))
                        quantum_state =(countj, countka, countkc, countm, args.isomer)
                        # generate acceleration file from Stark-effect file and quantum state, and deflection field norm and gradient
                        acceleration_fields = generate_acceleration_field(quantum_state, stark_filename,
                                                                   deflector_fieldnorm_filename, deflector_fieldgradient_filename,
                                                                   deflector_voltage_scaling)
                        deflection_field, deflection_gradient_x, deflection_gradient_y = acceleration_fields
                        mueff = stark_effect(quantum_state, stark_filename)
                        # hit count
                        hit = 0
                        # open output file
                        table = output.create_table("/"+("_"+str(countj)+"_"+str(countka)+"_"+str(countkc)+"_"+str(countm)+"_"+str(args.isomer)), 'FlightData', Particle)
                        # let them fly
                        for n in range(num_molecules):
                            initial = sample_source(source_position, source_width)
                            final, time = fly(initial, (0, deflector_length), (skimmer1_position, skimmer1_width), (skimmer2_position, skimmer2_width), (skimmer3_position, skimmer3_width),detector_position, 0.001)
                            if final[2] >= detector_position:
                                # trajectory reached end of deflector -> write initial and final positions to output file
                                # save final phasespace in output file
                                table.append([(str(time), final[0],final[1],final[2],final[3],final[4],final[5],str(0.),initial[0],initial[1],initial[2],initial[3],initial[4],initial[5])])
                                table.flush()
                            else:
                                hit += 1
                        print("hit "+str(hit))
                        # all molecules flown, close output and finish
    output.close()

































