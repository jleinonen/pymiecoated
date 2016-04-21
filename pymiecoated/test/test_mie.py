#! /usr/bin/env python

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

import unittest
from ..mie_coated import Mie
import sys


#some allowance for rounding errors etc
epsilon = 1e3*sys.float_info.epsilon


def run_tests():
    """Tests for the Mie code.

       Runs several tests that test the Mie code. All tests should return ok.
       If they don't, please contact the author.
    """
    suite = unittest.TestLoader().loadTestsFromTestCase(MieTests)
    unittest.TextTestRunner(verbosity=2).run(suite)


class MieTests(unittest.TestCase):

    def test_single_nonmagnetic(self):
        mie = Mie(m=complex(1.5,0.5),x=2.5)

        qext_ref = 2.562873497454734
        qsca_ref = 1.0970718190883924
        qabs_ref = 1.4658016783663417
        qb_ref = 0.12358646817981821
        asy_ref = 0.74890597894850719
        qratio_ref = 0.112651210275834

        for (func,ref) in zip(
            (mie.qext, mie.qsca, mie.qabs, mie.qb, mie.asy, mie.qratio),
            (qext_ref,qsca_ref,qabs_ref,qb_ref,asy_ref,qratio_ref)):
            self.assertLess(abs(ref-func())/ref, epsilon)

        S12_ref = (complex(-0.49958438416709694,-0.24032581667666403),
                   complex(0.11666852712178288,0.051661382367147853))
        S12 = mie.S12(-0.6)
        self.assertLess(abs(S12_ref[0]-S12[0])/S12_ref[0], epsilon)
        self.assertLess(abs(S12_ref[1]-S12[1])/S12_ref[1], epsilon)


    def test_coated(self):
        mie = Mie(m=complex(1.5,0.5),m2=complex(1.2,0.2),x=1.5,y=5.0)

        qext_ref = 2.0765452928100769
        qsca_ref = 0.90777572021757091
        qabs_ref = 1.168769572592506
        qb_ref = 0.022692436240597712
        asy_ref = 0.90560220988567752
        qratio_ref = 0.024997844440209204

        for (func,ref) in zip(
            (mie.qext, mie.qsca, mie.qabs, mie.qb, mie.asy, mie.qratio),
            (qext_ref,qsca_ref,qabs_ref,qb_ref,asy_ref,qratio_ref)):
            self.assertLess(abs(ref-func())/ref, epsilon)

        S12_ref = (complex(0.28677219451960079,-0.063605895700765691),
                   complex(-0.32635924647084191,0.12670342074119806))
        S12 = mie.S12(-0.6)
        self.assertLess(abs(S12_ref[0]-S12[0])/S12_ref[0], epsilon)
        self.assertLess(abs(S12_ref[1]-S12[1])/S12_ref[1], epsilon)


    def test_magnetic(self):
        mie = Mie(eps=complex(2.2,0.8),mu=complex(1.6,1.4),x=4.0)

        qext_ref = 2.6665582594291073
        qsca_ref = 1.1255460946883893
        qabs_ref = 1.541012164740718
        qb_ref = 0.0072453174040961301
        asy_ref = 0.89955981937838192
        qratio_ref = 0.0064371574281033928

        for (func,ref) in zip(
            (mie.qext, mie.qsca, mie.qabs, mie.qb, mie.asy, mie.qratio),
            (qext_ref,qsca_ref,qabs_ref,qb_ref,asy_ref,qratio_ref)):
                
            self.assertLess(abs(ref-func())/ref, epsilon)

        S12_ref = (complex(0.14683196954000932,-0.017479181764394575),
                   complex(-0.12475414168001844,0.28120475717321358))
        S12 = mie.S12(-0.6)
        self.assertLess(abs(S12_ref[0]-S12[0])/S12_ref[0], epsilon)
        self.assertLess(abs(S12_ref[1]-S12[1])/S12_ref[1], epsilon)


    def test_keywords(self):
        mie = Mie(m=complex(1.5,0.5),m2=complex(1.2,0.2),x=1.5,y=5.0)
        mie2 = Mie()
        mie2.m = complex(1.5,0.5)
        mie2.m2 = complex(1.2,0.2)
        mie2.x = 1.5
        mie2.y = 5.0

        for (func1,func2) in zip(
            (mie.qext, mie.qsca, mie.qabs, mie.qb, mie.asy, mie.qratio),
            (mie2.qext, mie2.qsca, mie2.qabs, mie2.qb, mie2.asy, 
            mie2.qratio)):
            
            self.assertEqual(func1(),func2())


    def test_errors(self):
        mie = Mie()

        #test that negative values of x fail
        def test_x():
            mie.x = -1.0
        self.assertRaises(ValueError, test_x)
        mie.x = 1.0

        #test that a missing m fails
        self.assertRaises(ValueError, mie.qext)
        mie.m = complex(1.5,0.5)

        #test that y<x fails (y==x is permitted)
        def test_y():
            mie.y = 0.5
        self.assertRaises(ValueError, test_y)
        mie.y = 1.5

        #test that setting y without m2 fails
        self.assertRaises(ValueError, mie.qext)
        mie.m2 = complex(1.2,0.5)

        #test that invalid values of u fail
        self.assertRaises(ValueError, mie.S12, -1.5)

        mie.mu = complex(1.5,0.6)
        #test that multilayered particles with mu fail
        self.assertRaises(ValueError, mie.qext)


if __name__ == '__main__':
    unittest.main()
