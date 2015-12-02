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


class Detector:
    """Default detector implementation

    This detector takes the (final) phase position passed and saves it to the specified single HDF5 file.

    ..todo:: Should close the file upon destruction of the object
    """

    def __init__(self, outname, z_max, t_max):
        # create a pytables HDF5 file to store all final phase-space points in
        # should enforce prefix ".fly" (simply add if not there)
        self.__z_max = z_max
        self.__t_max = t_max


    def append(self, position):
        # add position to HDF5 table file
        pass


    def z_max(self): return self.__z_max


    def t_max(self): return self.__t_max
