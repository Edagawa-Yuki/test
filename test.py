#! c:/Python27/python.exe
# -*- coding: utf-8 -*-

import cmath
import numpy
from pylab import *
from cmath import sqrt


#xi = 33.1191585
#xi = 147.081297

start = 20
end = 35
period = 0.05
s=[]
for xi in arange(start, end, period):
    zeta1 = (1.0+3.0**(1.0/2.0)*1j) * (-3*xi -1) /3/4 ** (1.0/3.0) / (3.0*3.0**(1.0/2.0)*sqrt(-4.0*xi**3+71*xi**2-8*xi)-45*xi+2.0)**(1.0/3.0) - (1.0 - 3.0 **(1.0/2.0) * 1j) * (3.0*3.0**(1.0/2.0)*sqrt(-4*xi**3+71*xi**2-8*xi)-45*xi+2.0)**(1.0/3.0)/6.0/2.0**(1.0/3.0)-2.0/3.0
    zeta2 = (1.0-3.0**(1.0/2.0)*1j) * (-3*xi -1) /3/4 ** (1.0/3.0) / (3.0*3.0**(1.0/2.0)*sqrt(-4.0*xi**3+71*xi**2-8*xi)-45*xi+2.0)**(1.0/3.0) - (1.0 + 3.0 **(1.0/2.0) * 1j) * (3.0*3.0**(1.0/2.0)*sqrt(-4*xi**3+71*xi**2-8*xi)-45*xi+2.0)**(1.0/3.0)/6.0/2.0**(1.0/3.0)-2.0/3.0
    zeta3 = -2.0**(1.07/3.0)*(-3.0*xi-1.0)/(3.0*(3.0*sqrt(3.0)*sqrt(-4.0*xi**3+71*xi**2-8.0*xi)-45*xi+2.0)**(1.0/3.0))+(3.0*sqrt(3.0)*sqrt(-4.0*xi**3+71*xi**2-8.0*xi)-45*xi+2.0)**(1.0/3.0)/3.0/2.0**(1.0/3.0)-2.0/3.0
    d = sqrt(-1.0 * zeta1.real)
    e = sqrt(zeta2.real)
    f = sqrt(zeta3.real)
    #zeta = 0.828
    #mu = 0.509
    B = numpy.matrix([[1.0,0.0,1.0,0.0,1.0,0.0],[0.0,d,0.0,e,0.0,f],[-d**2,0.0,e**2,0.0,f**2,0.0],[cos(d*pi),sin(d*pi),cosh(e*pi),sinh(e*pi),cosh(f*pi),sinh(f*pi)],[-d**2*cos(d*pi),-d**2*sin(d*pi),e**2*cosh(e*pi),e**2*sinh(e*pi),f**2*cosh(f*pi),f**2*sinh(f*pi)],[d**4*cos(d*pi),d**4*sin(d*pi),e**4*cosh(e*pi),e**4*sinh(e*pi),f**4*cosh(f*pi),f**4*sinh(f*pi)]])
    s.append(numpy.linalg.det(B))
#A = numpy.matrix([[1.0,0.0,1.0,0.0,0.0,0.0],[0.0,d,0.0,mu,0.0,zeta],[-d**2,0.0,zeta**2-mu**2,0.0,2.0*zeta*mu,0.0],[cos(d*pi),sin(d*pi),cos(mu*pi)*cosh(zeta*pi),sin(mu*pi)*cosh(zeta*pi),cos(mu*pi)*sinh(zeta*pi),sin(mu*pi)*sinh(zeta*pi)],[-d**2*cos(d*pi),-d**2*sin(d*pi),(zeta**2-mu**2)*cos(mu*pi)*cosh(zeta*pi)-2*zeta*mu*sin(mu*pi)*sinh(zeta*pi),(zeta**2-mu**2)*sin(mu*pi)*cosh(zeta*pi)+2*zeta*mu*cos(mu*pi)*sinh(zeta*pi),(zeta**2-mu**2)*cos(mu*pi)*sinh(zeta*pi)+2*zeta*mu*sin(mu*pi)*cosh(zeta*pi),(zeta**2-mu**2)*sin(mu*pi)*sinh(zeta*pi)+2*zeta*mu*cos(mu*pi)*cosh(zeta*pi)],[d**4*cos(d*pi),d**4*sin(d*pi),(mu**4-6*zeta**2*mu**2+zeta**4)*cos(mu*pi)*cosh(zeta*pi)+4*mu*zeta*(mu**2-zeta**2)*sin(mu*pi)*sinh(zeta*pi),(mu**4-6*zeta**2*mu**2+zeta**4)*sin(mu*pi)*cosh(zeta*pi)-4*mu*zeta*(mu**2-zeta**2)*cos(mu*pi)*sinh(zeta*pi),(mu**4-6*zeta**2*mu**2+zeta**4)*cos(mu*pi)*sinh(zeta*pi)+4*mu*zeta*(mu**2-zeta**2)*sin(mu*pi)*cosh(zeta*pi),(mu**4-6*zeta**2*mu**2+zeta**4)*sin(mu*pi)*sinh(zeta*pi)-4*mu*zeta*(mu**2-zeta**2)*cos(mu*pi)*cosh(zeta*pi)]])
#A = numpy.matrix([[(mu**4-6*zeta**2*mu**2+zeta**4)*cos(mu*pi)*sinh(zeta*pi)+4*mu*zeta*(mu**2-zeta**2)*sin(mu*pi)*cosh(zeta*pi),(mu**4-6*zeta**2*mu**2+zeta**4)*sin(mu*pi)*sinh(zeta*pi)-4*mu*zeta*(mu**2-zeta**2)*cos(mu*pi)*cosh(zeta*pi)]])
#print B
t=arange(start,end,period)
plot(t,s)
show()
