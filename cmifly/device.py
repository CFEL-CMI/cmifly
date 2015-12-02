# -*- coding: utf-8; mode: python; fill-column: 120; truncate-lines: t -*-
# cython: language_level=3
#
# This file is part of CMIfly
# Copyright (C) 2015 Jochen Küpper <jochen.kuepper@cfel.de>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this programm for scientific work, you must correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

__author__ = "Jochen Küpper <jochen.kuepper@cfel.de>"

import numpy as np
import tables


class Device:
    """Device description"""
    pass


class Deflector_A(Device):
    """The electric deflector (a-type geometry)

    The a-type deflector is constructed as a 3D object that has a transverse mechanical aperture between two circles, a
    2D transverse electric field norm and gradient (stored in a PyTables object), and an extend from z = z1 to z2.
    """

    def __init__(self, filename, z=(0.,0.1)):
        """Construct the deflector"""
        # should actually
        self.__field = tables.open_file(filename, mode='r')
        self.__z = np.array(z)
        assert((2,) == self.__z.shape)



class Deflector_B(Device):
    """The electric deflector (b-type geometry)

    The b-type deflector is constructed as a 3D object that has a transverse mechanical aperture of a square, a
    2D transverse electric field norm and gradient (stored in a PyTables object), and an extend from z = z1 to z2.
    """

    def __init__(self, filename, z=(0.,0.1)):
        """Construct the deflector"""
        # should actually
        self.__field = tables.open_file(filename, mode='r')
        self.__z = np.array(z)
        assert((2,) == self.__z.shape)



class LaserBeam(Device):
    """Description of the intensity and polarization distribution of a laser beam in space"""
    pass



class GasPhase_environment(Device):
    """Describe the properties of a gas-phase environment

    Currently, all gas particles are the same and the pressure and the temperature is the same across all space
    """

    def __init__(self, pressure, temperature, particle):
        self.__p = pressure
        self.__T = temperature
        self.__particle = particle


    def pressure(pos):
        """Pressure at position pos

        pos is a 3-tuple or np.array of shape (3,)
        """
        return self.__p

    def temperature(pos):
        """Temperature at position pos

        pos is a 3-tuple or np.array of shape (3,)
        """
        return self.__T

    def particle_mass(self):
        return self.__particle.mass()



class Laser_Gas_Interaction(Device):
    """Laser-gas interaction

    Description of the inetraction of a particle with a laser field in a gaseous medium.

    This model takes into account the optical and thermal properties of a Particle, the photophoretic force of the laser
    field onto that particle in the given gas environment, and the

    """

    def __init__(self, laserbeam, environment):
        """Construct the device"""
        # should actually
        self.__field = tables.open_file(filename, mode='r')
        self.__z = np.array(z)
        assert((2,) == self.__z.shape)
