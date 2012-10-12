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

from numpy import pi, arange, zeros, hstack, sqrt, sin, cos
from scipy.special import jv, yv


class MieCoeffs(object):   
   """Wrapper for the Mie coefficients.
   """
   def __init__(self, par):         
      (self.an, self.bn, self.nmax) = mie_coeffs(par)
      

def mie_coeffs(params):
   """Input validation and function selection for the Mie coefficients.
   """
   
   eps = complex(params["eps"]) if params["eps"] is not None else None
   x = float(params["x"]) if params["x"] is not None else None
   
   if (x==None) or (eps==None):
      raise ValueError("Must specify x and either eps or m.")
   
   mu = complex(params["mu"]) if params["mu"] is not None else complex(1.0)
   
   y = float(params["y"]) if params["y"] is not None else None
   eps2 = complex(params["eps2"]) if params["eps2"] is not None else None 
         
   coated = (y is not None)  
   if coated == (eps2 is None):
      raise ValueError("Must specify both y and m2 for coated particles.")
      
   if coated and mu!=complex(1.0):
      raise ValueError("Multilayer calculations for magnetic particles are not currently supported.")
   
   if not coated:
      y = x
      eps2 = eps   
     
   # Do not use the coated version if it is not necessary   
   if x==y or eps==eps2:      
      coeffs = single_mie_coeff(eps,mu,y)
   elif x==0:   
      coeffs = single_mie_coeff(eps2,mu,y)
   else:
      coeffs = coated_mie_coeff(eps,eps2,x,y)
      
   return coeffs


def single_mie_coeff(eps,mu,x):
   """Mie coefficients for the single-layered sphere.
   
      Args:
         eps: The complex relative permittivity.
         mu: The complex relative permeability.
         x: The size parameter.
         
      Returns:
         A tuple containing (an, bn, nmax) where an and bn are the Mie
         coefficients and nmax is the maximum number of coefficients.         
   """
   z = sqrt(eps*mu)*x
   m = sqrt(eps/mu)

   nmax = int(round(2+x+4*x**(1.0/3.0)))
   nmax1 = nmax-1
   nmx = int(round(max(nmax,abs(z))+16))
   n = arange(nmax)
   nu = n+1.5
   
   sx = sqrt(0.5*pi*x)
   px = sx*jv(nu,x)
   p1x = hstack((sin(x), px[:nmax1]))
   chx = -sx*yv(nu,x)
   ch1x = hstack((cos(x), chx[:nmax1]))
   gsx = px-complex(0,1)*chx
   gs1x = p1x-complex(0,1)*ch1x
   
   dnx = zeros(nmx,dtype=complex)
   for j in xrange(nmx-1,0,-1):
      r = (j+1.0)/z
      dnx[j-1] = r - 1.0/(dnx[j]+r)
   dn = dnx[:nmax]
   n1 = n+1
   da = dn/m + n1/x
   db = dn*m + n1/x
   
   an = (da*px-p1x)/(da*gsx-gs1x)
   bn = (db*px-p1x)/(db*gsx-gs1x)
   
   return (an, bn, nmax)   


def coated_mie_coeff(eps1,eps2,x,y):
   """Mie coefficients for the dual-layered (coated) sphere.
   
      Args:
         eps: The complex relative permittivity of the core.
         eps2: The complex relative permittivity of the shell.
         x: The size parameter of the core.
         y: The size parameter of the shell.
         
      Returns:
         A tuple containing (an, bn, nmax) where an and bn are the Mie
         coefficients and nmax is the maximum number of coefficients.         
   """   
   m1 = sqrt(eps1)
   m2 = sqrt(eps2)
   m = m2/m1
   u = m1*x
   v = m2*x
   w = m2*y
   
   nmax = int(round(2+y+4*y**(1.0/3.0)))
   mx = max(abs(m1*y),abs(w))
   nmx = int(round(max(nmax,mx)+16))
   nmax1 = nmax-1
   n = arange(nmax)
   
   dnu = zeros(nmax,dtype=complex)
   dnv = zeros(nmax,dtype=complex)
   dnw = zeros(nmax,dtype=complex)
   dnx = zeros(nmx,dtype=complex)   
   
   for (z, dn) in zip((u,v,w),(dnu,dnv,dnw)):      
      for j in xrange(nmx-1,0,-1):
         r = (j+1.0)/z
         dnx[j-1] = r - 1.0/(dnx[j]+r)
      dn[:] = dnx[:nmax]
   
   nu = n+1.5
   vwy = [v,w,y]
   sx = [sqrt(0.5*pi*xx) for xx in vwy]
   (pv,pw,py) = [s*jv(nu,xx) for (s,xx) in zip(sx,vwy)]
   (chv,chw,chy) = [-s*yv(nu,xx) for (s,xx) in zip(sx,vwy)]   
   p1y = hstack((sin(y), py[:nmax1]))
   ch1y = hstack((cos(y), chy[:nmax1]))
   gsy = py-complex(0,1)*chy
   gs1y = p1y-complex(0,1)*ch1y
   
   uu = m*dnu-dnv
   vv = dnu/m-dnv   
   fv = pv/chv
   fw = pw/chw
   ku1 = uu*fv/pw
   kv1 = vv*fv/pw
   pt = pw-chw*fv
   prat = pw/pv/chv
   ku2 = uu*pt+prat
   kv2 = vv*pt+prat
   dns1 = ku1/ku2
   gns1 = kv1/kv2

   dns = dns1+dnw
   gns = gns1+dnw
   nrat = (n+1)/y
   a1 = dns/m2+nrat
   b1 = m2*gns+nrat
   an = (py*a1-p1y)/(gsy*a1-gs1y)
   bn = (py*b1-p1y)/(gsy*b1-gs1y)
   
   return (an, bn, nmax)
