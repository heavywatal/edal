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
F_approx = beta * triangle

Xi_approx = 1 - h0 * ((u - y0) ** 2) - h1 * ((v - y1) ** 2)
Xi_normal = exp( - h0 * ((u - y0) ** 2) - h1 * ((v - y1) ** 2))
T = F_approx * Xi_approx


with assuming(Q.positive(a - 1)):
    print('Xi_approx triangle')
    xi_denom = integrate(integrate(Xi_approx, (v, 0, 1 - u)), (u, 0, 1))
    print(xi_denom)
    print(xi_denom.simplify())
    print(xi_denom.factor())

    print('Xi_approx square')
    xi_denom = integrate(integrate(Xi_approx, (v, 0, 1)), (u, 0, 1))
    print(xi_denom)
    print(xi_denom.simplify())
    print(xi_denom.factor())

    if False:
        print('Xi_normal')
        print(integrate(Xi_normal, (v, 0, 1)).simplify())

        print('D_I')
        D_I = integrate(integrate(T, (v, 0, 1 - u)), (u, 0, 1))
        print(D_I.simplify())

