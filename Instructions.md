# Installation #
Make sure you have [NumPy](http://numpy.org/) and [SciPy](http://scipy.org/) (and of course, a Python environment) installed.

The code is easiest to install from PyPI. If you have  [pip](http://pypi.python.org/pypi/pip) installed, install simply with
```
$ pip install pymiecoated
```
(using sudo if installing globally).

If you can't use `pip`, install manually:
  1. [Download](http://code.google.com/p/pymiecoated/downloads/list) the latest package and unzip OR [clone](http://code.google.com/p/pymiecoated/source/checkout) the latest version from the repository.
  1. Go to the root directory of the package and run (with sudo if you need the permissions):
```
$ python setup.py install 
```


# Testing #

This is to verify that you have a working package after installation. Open a Python shell. Run the following:
```
from pymiecoated.test import test_mie
test_mie.run_tests()
```
An output like this should be printed:
```
test_coated (pymiecoated.test.test_mie.MieTests) ... ok
test_errors (pymiecoated.test.test_mie.MieTests) ... ok
test_keywords (pymiecoated.test.test_mie.MieTests) ... ok
test_magnetic (pymiecoated.test.test_mie.MieTests) ... ok
test_single_nonmagnetic (pymiecoated.test.test_mie.MieTests) ... ok

----------------------------------------------------------------------
Ran 5 tests in 0.017s

OK
```
The running time may vary; the most important thing is that "ok" is printed after each test. If some of the tests fail, please notify the author.

Please note that these tests only verify the basic functionality of the package, they do not (as of the latest version) test for convergence in edge cases.

# Usage #

## Minimal example ##
This calculates the scattering efficiency for a sphere with size parameter 1.5 and refractive index 3.0+0.5i:
```
>>> from pymiecoated import Mie
>>> mie = Mie(x=1.5,m=complex(3.0,0.5))
>>> mie.qsca()
1.4319737923076081
```

## Attributes ##
The Mie class object has a number of attributes that define the scatterer. They can be set either using keyword arguments when creating the object, or by setting object attributes directly. This code
```
>>> mie = Mie(x=1.5,m=complex(3.0,0.5))
```
should be equal to this:
```
>>> mie = Mie()
>>> mie.x = 1.5
>>> mie.m = complex(3.0,0.5)
```

The following attributes are recognized:
| _x_ | The size parameter of the sphere, or 2\*pi\*r/lambda for sphere radius r and wavelength lambda. If _y_ is defined, this is the size parameter of the inner layer (core).|
|:----|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _eps_ | The complex relative permittivity (dielectric constant). If _eps2_ is defined, this is the relative permittivity of the inner layer.                                    |
| _m_ | The complex refractive index. If m2 is defined, this is the refractive index of the core.                                                                               |
| _mu_ | The complex relative permeability.                                                                                                                                      |
| _y_ | The size parameter of the outer layer (shell).                                                                                                                          |
| _eps2_ | The complex relative permittivity of the outer layer.                                                                                                                   |
| _m2_ | The complex refractive index of the outer layer.                                                                                                                        |

The attributes _m_ and _m2_ are special in the sense that they are shortcuts to the attributes _eps_ and _eps2_. Setting _m_ will also set _mu_ to zero. Currently, multilayered magnetic spheres are not supported, so setting _mu_ together with _m2_, _eps2_ or _y_ will raise an error.

## Methods ##
These methods return the scattering parameters:
|_qext()_|Extinction efficiency (get the cross section by multiplying the efficiency with pi\*r^2 where r is the radius of the sphere)|
|:-------|:---------------------------------------------------------------------------------------------------------------------------|
|_qsca()_|Scattering efficiency                                                                                                       |
|_qabs()_|Absorption efficiency                                                                                                       |
|_qb()_  |Backscattering efficiency                                                                                                   |
|_asy()_ |Asymmetry parameter                                                                                                         |
|_qratio()_|The backscattering ratio, i.e. _qb()_/_qsca()_.                                                                             |
|_S12(u)_|The amplitude scattering matrix according to the definition of Bohren and Huffman (1983). The parameter _u_ (-1<=_u_<=1) is the cosine of the scattering angle.|

## Example (magnetic sphere) ##
A magnetic sphere with a size parameter 1.5, relative permittivity 5.5+0.8i, and relative permeability 1.1+2.0i.
```
>>> mie = Mie(x=1.5,eps=complex(5.5,0.8),mu=complex(1.1,2.0))
>>> mie.qsca()
1.1197778395829314
```

## Example (coated sphere) ##
A dual-layered sphere with a core size parameter 1.5, core refractive index 3.0+0.5i, shell size parameter 3.0, and shall refractive index 2.0+0.2i.
```
>>> mie = Mie(x=1.5,m=complex(3.0,0.5),y=3.0,m2=complex(2.0,0.2))
>>> mie.S12(-0.6)
((-0.16412698606840209+0.81678948249170791j),
 (1.0709433632834138-0.12417452542839993j))
```

# Code documentation #
The public interface of the code is well documented with Python docstrings (the internals, less so). If you have [IPython](http://ipython.org/), you can type a question mark after an object to get the docstring. For example:
```
In [1]: from pymiecoated import Mie
In [2]: mie = Mie()
In [3]: mie.qb?
Type:       instancemethod
Base Class: <type 'instancemethod'>
String Form:<bound method Mie.qb of <pymiecoated.mie_coated.Mie object at 0x3cd2e90>>
Namespace:  Interactive
Definition: mie.qb(self)
Docstring:
Returns:
The backscattering efficiency. Multiply by the physical cross section
(pi*r**2) to get the cross section.
```

# Demo #
A demo that computes backscattering cross sections for spherical hailstones at various stages of melting is can be run as follows:
```
>>> from pymiecoated.demos import mie_demos
>>> mie_demos.melting_hail()
```
This demo requires [matplotlib](http://matplotlib.org/) to run.

# Citing the code #
No paper publication to cite is available. If you'd like to cite this code in a publication, the following reference is suggested:

Leinonen, J., _Python code for calculating Mie scattering from single- and dual-layered spheres_. Available at http://code.google.com/p/pymiecoated/.

# Original code #
The original MATLAB code by C. Mätzler is available at: http://diogenes.iwt.uni-bremen.de/vt/laser/codes/Mie-Matlab-Maetzler.zip

# References #
Bohren, C. F. and D. R. Huffman (1983), _Absorption and Scattering of Light by Small Particles_. John Wiley & Sons.