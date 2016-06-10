"""
Copyright (C) 2012-2013 Jussi Leinonen

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

from numpy import sqrt
from mie_coeffs import MieCoeffs
from mie_aux import Cache
from mie_props import mie_props, mie_S12


class MieScatterProps(object):
    """Stores the mie coefficients and the corresponding parameters.
    """
    def __init__(self, params):
        par = dict(zip(("eps","mu","x","y","eps2"),params[:5]))
        if par["x"]==0 and par["y"] is None:
            #give valid output for x==0
            self._coeffs = {"qext":0.0, "qsca":0.0, "qabs":0.0, "qb":0.0,
                            "0.0":asy, "0.0":qratio}
        else:
            self._coeffs = MieCoeffs(par)
        self._props = None
        self._S12 = None
        self.size = par["x"] if par["y"]==None else par["y"]

    def prop(self, prop_name):
        if self._props is None:
            self._props = mie_props(self._coeffs, self.size)
        return self._props[prop_name]

    def S12(self, u):
        self._S12 = mie_S12(self._coeffs, u)
        return self._S12


class Mie(object):
    """Class for computing Mie scattering from homogeneous and coated spheres.

    Create an instance of Mie, set the size and material parameters, then call
    any of the methods to determine the scattering properties.

    Attributes:
    x: The size parameter, i.e. 2*pi/lambda*r where lambda is the
        wavelength and r is the radius of the spheres.
    eps: The complex relative permittivity (dielectric constant).
    mu: The complex relative permeability.
    eps2: The dielectric constant of the outer layer. If you specify this,
        you must also specify y. If you specify eps2, eps is taken to be
        the complex relative permittivity of the inner layer (core).
    y: The size parameter of the outer layer. If you specify y, x is
        taken to be the size paramater of the core (must be y>=x).
    m: The complex refractive index. Setting this sets eps to m**2 and mu
        to 1.0.
    m2: The complex refractive index of the outer layer. See "m" above.

    Setting mu together with eps2 and y raises an error.

    Any of the above attributes can be given as keyword arguments when
    creating a new Mie instance. For example:
    mie = Mie(x=1.5,m=complex(1.2,0.1))
    """
    def __init__(self, **kwargs):
        self._cache = Cache()
        self.eps = None
        self.mu = 1.0
        self._x = None
        self._y = None
        self.eps2 = None
        for k in ["eps","mu","eps2"]:
            if k in kwargs:
                self.__dict__[k] = kwargs[k]
        if "m" in kwargs:
            self.m = kwargs["m"]
        if "m2" in kwargs:
            self.m2 = kwargs["m2"]
        if "x" in kwargs:
            self.x = kwargs["x"]
        if "y" in kwargs:
            self.y = kwargs["y"]

    def _params_signature(self):
        return (self.eps, self.mu, self.x, self.y, self.eps2)

    def qext(self):
        """The extinction efficiency.

        Returns:
            The extinction efficiency. Multiply by the physical cross section
            (pi*r**2) to get the cross section.
        """
        return self._get_scatt_prop("qext")

    def qsca(self):
        """The scattering efficiency.

        Returns:
            The scattering efficiency. Multiply by the physical cross section
            (pi*r**2) to get the cross section.
        """
        return self._get_scatt_prop("qsca")

    def qabs(self):
        """The absorption efficiency.

        Returns:
            The absorption efficiency. Multiply by the physical cross section
            (pi*r**2) to get the cross section.
        """
        return self._get_scatt_prop("qabs")

    def qb(self):
        """The backscattering efficiency.

        Returns:
            The backscattering efficiency. Multiply by the physical cross
            section (pi*r**2) to get the cross section.
        """
        return self._get_scatt_prop("qb")

    def asy(self):
        """The asymmetry parameter.
        Returns:
            The asymmetry parameter, i.e. <cos(theta)>.
        """
        return self._get_scatt_prop("asy")

    def qratio(self):
        """The backscattering ratio.

        Returns:
              The backscattering ratio, i.e. qb()/qsca().
        """
        return self._get_scatt_prop("qratio")

    def S12(self, u):
        """The amplitude scattering matrix elements.

        Arguments:
            u: The cosine of the scattering angle, -1 <= u <= 1

        Returns:
            The amplitude scattering matrix elements S1 and S2.
            Follows the conventions of Bohren and Huffman (1983).
        """
        return self._get_S12(u)


    def _get_scatt_prop(self, prop):
        sig = self._params_signature()
        if sig not in self._cache:
            self._cache[sig] = MieScatterProps(sig)
        return self._cache[sig].prop(prop)

    def _get_S12(self, u):
        if abs(u) > 1:
            raise ValueError("The cosine u must be between -1 and 1.")
        sig = self._params_signature()
        if sig not in self._cache:
            self._cache[sig] = MieScatterProps(sig)
        return self._cache[sig].S12(u)


    def _get_m(self):
        return sqrt(self.eps/self.mu)

    def _set_m(self, m):
        self.mu = 1.0
        self.eps = m**2

    m = property(_get_m, _set_m)


    def _get_m2(self):
        return sqrt(self.m2)

    def _set_m2(self, m2):
        self.eps2 = m2**2

    m2 = property(_get_m2, _set_m2)


    def _get_x(self):
        return self._x

    def _set_x(self, x):
        if x >= 0.0:
            self._x = x
        else:
            raise ValueError("The size x cannot be smaller than 0.")

    x = property(_get_x, _set_x)


    def _get_y(self):
        return self._y

    def _set_y(self, y):
        if y >= self.x:
            self._y = y
        else:
            raise ValueError("The size y cannot be smaller than x.")

    y = property(_get_y, _set_y)
