#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2008,2009,2012,2013,2014,2015 Jochen Küpper <jochen.kuepper@cfel.de>


import numpy,os
from setuptools import setup
from Cython.Build import cythonize

extra_compile_args = []
library_dirs = []

long_description = """CMIfly -- trajectory calculations for individual particles in force fields

Developed by the Controlled Molecule Imaging Group at the Center for Free-Electron Laser Science, DESY and Universität
Hamburg, Germany.

Original author:     Jochen Küpper <jochen.kuepper@cfel.de>
Current maintainer:  Jochen Küpper <jochen.kuepper@cfel.de>
See the distribution files AUTHORS and THANKS for further contributions.
"""


setup(name="cmifly",
      author              = "Jochen Küpper, CFEL-CMI group, et al (see AUTHORS)",
      author_email        = "jochen.kuepper@cfel.de",
      maintainer          = "Jochen Küpper and the CFEL-CMI group",
      maintainer_email    = "jochen.kuepper@cfel.de",
      url                 = "https://www.controlled-molecule-imaging.org",
      description         = "CMIfly -- trajectory calculations for individual particles in force fields",
      version             = "0.1.dev0",
      long_description    = long_description,
      license             = "GPL",
      packages            = ['cmifly'],
      scripts             = ['scripts/cmifly',
                             'scripts/cmifly_plot-phasespace',
                             'scripts/cmifly_thermally-averaged-deflection'
      ],
      test_suite          = 'tests',
      )
