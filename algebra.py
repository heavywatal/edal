#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
"""
from sympy import *

u, v = symbols('u v', real=True, positive=True)
y0, y1 = symbols('y0 y1', real=True, positive=True)
h0, h1 = symbols('h0 h1', real=True, positive=True)
a = Symbol('a', real=True, positive=True)
f = 2 / (1 - u) * (1 - (u / (1 - v)))
g = (u * (1 - u)) ** a
h = 1 - h0 * (u - y0) - h1 * (v - y1)
F = f * g * h
print(F.simplify())

with assuming(Q.positive(a - 1)):
    G = integrate(F, (v, 0, 1 - u))
    print(G.simplify())

#    H = integrate(G, (u, 0, 1))
#    print(H)

#http://integrals.wolfram.com/index.jsp
#(2*x**(a + 1)*(-x + 1)**a*(-h0*x + h0*y0 + h1*y1 - h1 + 1)*Log(-x + 1) - 2*x**(a + 1)*(-x + 1)**a*(-h0*x + h0*y0 + h1*y1 - h1 + 1)*Log(-x - (x - 1)**2 + 1) + (x*(-x + 1))**a*(x - 1)*(-2*h0*x + 2*h0*y0 - 2*h1*x + 2*h1*y1 + h1*(x - 1) + 2))/(x - 1)
