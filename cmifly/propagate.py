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

__doc__ = """Propagate particles in force field"""

import numpy as np


def propagate(particle, accel, detector):
    """particle trajectory calculation integrator

    Implement a particle trajectory calculation integrator that transports the particle from z0/t0 through the
    acceleration field until one of the position z_max or time t_max from the detector are reached.

    The particles' initial positions include z_0 and t_0

    Returns the final phase-space position.

    """`
    pass


def fly(source, acceleration, detector, n=1e6):
    """Calculate trajectories for n particles"""
    for range(n):
        propagate(source.next(), acceleration, detector)
