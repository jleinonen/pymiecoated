#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

long_description = """A Python code for computing the scattering properties
of single- and dual-layered spheres with an easy-to-use object oriented
interface.

Based on code by C. Mätzler; ported and published with permission.

Requires NumPy and SciPy.
"""

setup(name='pymiecoated',
      version='0.1.2',
      download_url=\
          'http://pymiecoated.googlecode.com/files/pymiecoated-0.1.2.zip',
      description='Single- and dual-layered Mie scattering computations',
      author='Jussi Leinonen',
      author_email='jsleinonen@gmail.com',
      url='http://code.google.com/p/pymiecoated/',
      packages=['pymiecoated','pymiecoated.demos','pymiecoated.test'],
      license='MIT',
      long_description = long_description,
     )
