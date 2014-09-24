#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
http://integrals.wolfram.com/index.jsp
"""
from sympy import *

u, v = symbols('u v', real=True, positive=True)
y0, y1 = symbols('y0 y1', real=True, positive=True)
h0, h1 = symbols('h0 h1', real=True, positive=True)
a = Symbol('a', real=True, positive=True)

triangle = 2 / (1 - u) * (1 - (v / (1 - u)))
beta = (u * (1 - u)) ** (a - 1)
F = beta * triangle

Xi_approx = 1 - h0 * ((u - y0) ** 2) - h1 * ((v - y1) ** 2)
Xi_normal = exp( - h0 * ((u - y0) ** 2) - h1 * ((v - y1) ** 2))
T = F * Xi_approx


with assuming(Q.positive(a - 1)):
    print('D_I')
    D_I = integrate(T, (v, 0, 1 - u))
    print(D_I.simplify())

    print('Xi_approx')
    print(integrate(Xi_approx, (v, 0, 1 - u)).simplify())
    #(u - 1)*(3*h0*u**2 - 6*h0*u*y0 + 3*h0*y0**2 + 3*h1*y1**2 + 3*h1*y1*(u - 1) + h1*(u - 1)**2 - 3)/3

    print('Xi_normal')
    print(integrate(Xi_normal, (v, 0, 1 - u)).simplify())
