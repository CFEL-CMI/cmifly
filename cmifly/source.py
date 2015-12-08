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
import numpy.random

class Source:
    """Default Source implementation

    Provide a 6d phase-space position drawn from a uniform random number generator over [0,1]
    """

    def __init__(self):
        # we should pass in at least the center position and widths of the distributions
        pass


    def next(self, n=1):
        # return next $n$ molecule
        return np.random.uniform(size=6)

class CollimatedSource:
    """Collimated particle source

    Provide a 6d phase-space position. The x and y positionis drawn from a normal random number generator. z and the velocities are pre defined by the user.
    """

    def __init__(self, center, sigma, velocity ):
        self.__center = center
        self.__sigma = sigma
        self.__v = velocity

    def next(self)
        # return next molecule
        return [np.random.normal(self.__center[0],self.__sigma.[0],1),np.random.normal(self.__center[1],self.__sigma.[1],1),self.__center[2],0.,0.,-abs(self.__v)]

class FocusedSource:
    """Focused particle source
        
        Provide a 6d phase-space position. The x and y positionis drawn from a normal random number generator. z is pre defined by the user. The velocities are adapted to get a focus with user defined position and width.
        """
    
    def __init__(self, center, sigma, velocity, focus_z_position, focus_sigma ):
        self.__center = center
        self.__sigma = sigma
        self.__v = velocity
        self.__f = focus_z_position
        self.__sigma_f = focus_sigma
    
    
    def next(self)
    # return next molecule
    focus_distance=[np.random.normal(0,self.__sigma_f.[0],1),np.random.normal(0,self.__sigma_f.[1],1),-abs(self.__f)]
    t=np.sqrt(focus_distance[0]**2+focus_distance[1]**2+focus_distance[2]**2)/abs(self.__v)
    return [np.random.normal(self.__center[0],self.__sigma.[0],1),np.random.normal(self.__center[1],self.__sigma.[1],1),self.__center[2],focus_distance[0]/t,focus_distance[1]/t,focus_distance[2]/t]
