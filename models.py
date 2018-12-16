#!/usr/bin/env python3

from playground import analyze


def lorentz(a, value):
    # Example: ./models.py 0 .001 .999 1001 0 1e-9 1e-9 | ./plotMany.py 1 10 >/dev/null
    return (1 - a.sqr) ** -0.05 - 1 / (1 - a.sqr).sqrt - value


def cosx_x3(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a * a.sqr - value


def septic(a, value):
    # Example: ./models.py 2 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def composite1(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a.sqr - 4).exp).ln - value


def composite2(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
    return (a.sqr + (a.exp - 4).sqr).sqrt - value


def playground(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos - value
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr - value
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
    # return a.sqr.sqrt - value
    # return a.exp.ln - value
    # return a.tan - a.sin / a.cos - value
    # return a.tanh - a.sinh / a.cosh - value
    # return a.sin.asin - value
    # return a.cos.acos - value
    # return a.tan.atan - value
    # return a.asin - value
    # return a.acos - value
    return a.atan - value


analyze(composite1)
