#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Copyright (C) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import numpy
from numpy import array, pi
from matplotlib import pyplot
from ..mie_coated import Mie


def melting_hail():
   """
      A demo that calculates the backscattering cross sections of melting
      hail at different stages of melting. The hailstones are modeled with
      an ice core and a water shell.
   """
   c = 299792458.0
   wl = c/5.6e9*1e3
   m_i = complex(1.7844, 1.4773e-4)
   m_w = complex(8.3355, 2.2173)
   
   D_shell = numpy.linspace(1, 30, 1000)
   frac = numpy.arange(0.0,1.01,0.2)

   k = 2*pi/wl
   x_shell = k*D_shell
   
   xsect = {}
   
   mie = Mie(m=m_i,m2=m_w)
   def xsect_b(xc,xs):
         mie.x = xc
         mie.y = xs
         return mie.qb()*pi*xs**2
   
   for f in frac:      
      x_core = x_shell * f**(1.0/3.0)      
      xsect[f] = array([xsect_b(xc,xs) for (xc,xs) in zip(x_core,x_shell)])
      pyplot.semilogy(D_shell, xsect[f], label="Ice fraction "+str(f))
   
   pyplot.xlabel(u"Diameter (mm)", fontsize=14)
   pyplot.ylabel(u"Backscattering cross section (mm²)", fontsize=14)
   pyplot.legend(loc="lower right")
   
   return xsect
