#!/usr/bin/env python

## Program:   PyNS
## Module:    Bessel.py
## Language:  Python
## Date:      $Date: 2012/04/23 14:38:11 $
## Version:   $Revision: 0.4 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390


cdef extern from "complexobject.h":
 
    struct Py_complex:
        double real
        double imag
    
    ctypedef class __builtin__.complex [object PyComplexObject]:
        cdef Py_complex cval


def Bessel(int n, double complex arg):

    cdef double complex z, zproduct, zanswer, zarg
    cdef double k
    cdef int i
    
    z = 1. + 0.j
    zproduct = 1. + 0.j
    zanswer = 1. + 0.j
    zarg = -0.25 * (arg * arg)
    
    for i in range(0, 1000):
        k = (i+1.)*(i+1.+n)
        z = (1./(k))*(z*zarg)
        if abs(z) < 1e-20:
            break
        zanswer = zanswer + z
    for i in range(0,n):
        zproduct = zproduct * 0.5 * arg
    zanswer = zanswer * zproduct
    return zanswer