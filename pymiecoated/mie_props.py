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

from numpy import arange, dot, zeros, vstack, sqrt, sin, cos


def mie_props(coeffs,y):  
   """The scattering properties.
   """
   anp = coeffs.an.real
   anpp = coeffs.an.imag
   bnp = coeffs.bn.real
   bnpp = coeffs.bn.imag   
   nmax = coeffs.nmax
   
   n1 = nmax-1
   n = arange(1,nmax+1,dtype=float)
   cn = 2*n+1
   c1n = n*(n+2)/(n+1)
   c2n = cn/(n*(n+1))
   y2 = y**2
   
   dn = cn*(anp+bnp)
   q = dn.sum()
   qext = 2*q/y2
   
   en = cn*(anp**2+anpp**2+bnp**2+bnpp**2)
   q = en.sum()
   qsca = 2*q/y2
   qabs = qext-qsca   
   
   fn = (coeffs.an-coeffs.bn)*cn
   gn=(-1)**n
   q = (fn*gn).sum()
   qb = dot(q,q.conj()).real/y2
   
   g1 = zeros((4,nmax),dtype=float)
   g1[:,:n1] = vstack((anp[1:nmax], anpp[1:nmax], bnp[1:nmax], bnpp[1:nmax]))

   asy1 = c1n*(anp*g1[0,:]+anpp*g1[1,:]+bnp*g1[2,:]+bnpp*g1[3,:])
   asy2 = c2n*(anp*bnp+anpp*bnpp)
   
   asy = 4/y2*(asy1+asy2).sum()/qsca
   qratio = qb/qsca
   
   return {"qext":qext, "qsca":qsca, "qabs":qabs, "qb":qb, "asy":asy, "qratio":qratio}


def mie_S12(coeffs,u):  
   """The amplitude scattering matrix.
   """
   (pin,tin) = mie_pt(u,coeffs.nmax)
   n = arange(1, coeffs.nmax+1, dtype=float)
   n2 = (2*n+1)/(n*(n+1))
   pin *= n2
   tin *= n2   
   
   S1 = dot(coeffs.an,pin)+dot(coeffs.bn,tin)
   S2 = dot(coeffs.an,tin)+dot(coeffs.bn,pin)
   return (S1, S2)
   

def mie_pt(u,nmax):
   u = float(u)
   p = zeros(nmax, dtype=float)   
   p[0] = 1
   p[1] = 3*u
   t = zeros(nmax, dtype=float)
   t[0] = u
   t[1] = 6*u**2 - 3
   
   nn = arange(2,nmax,dtype=float)
   for n in nn:
      p[n] = (2*n+1)/n*p[n-1]*u - (n+1)/n*p[n-2]
      
   t[2:] = (nn+1)*u*p[2:] - (nn+2)*p[1:-1]
   
   return (p,t)
