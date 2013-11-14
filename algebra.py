#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
"""
from sympy import *

u, v = symbols('u v', real=True, positive=True)
y0, y1 = symbols('y0 y1', real=True, positive=True)
h0, h1 = symbols('h0 h1', real=True, positive=True)
a = Symbol('a', real=True, positive=True)
f = 2 / (1 - u) * (1 - (v / (1 - u)))
g = (u * (1 - u)) ** (a - 1)
h = 1 - h0 * ((u - y0) ** 2) - h1 * ((v - y1) ** 2)
F = f * g * h
print(F.simplify())

with assuming(Q.positive(a - 1)):
    G = integrate(F, (v, 0, 1 - u))
    print(G.simplify())

#    H = integrate(G, (u, 0, 1))
#    print(H)

#http://integrals.wolfram.com/index.jsp
