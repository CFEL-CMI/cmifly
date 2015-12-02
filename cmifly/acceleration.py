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

__doc__ = """Calculation of the acceleration field"""

import numpy as np
import tables


class Acceleration:
    """Implementation of Acceleration for the trivial case of field-free space"""

    def __init__(self):
        """Initialize the Acceleration field from the relevant parameters -- none in this field-free case"""
        pass

    def acceleration(x, y, z):
        """Calculate acceleration at position (x,y,z)"""
        return self.__getitem__((x,y,z))

    def acceleration(pos):
        """Calculate acceleration at position pos=(x,y,z)"""
        return self.__getitem__(pos)

    def __getitem__(self, pos):
        """Calculate acceleration at position pos=(x,y,z)"""
        return 0.;

    def __setitem__(self, key, value):
        pass



class DeflectorAcceleration(Acceleration):
    pass



class PhotoAcceleration(Acceleration):
    """Implementation of the acceleration of a "classical" particle in a laser beam.

    Includes the following forces
    * photophoretoc force according to (reference)
    * photon pressure

    Other effects, such as gvravity, buoyancy, viscosity, drag forces are neglected so far.
    """
    pass
