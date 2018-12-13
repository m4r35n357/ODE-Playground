#!/usr/bin/env python3

from playground import analyze


def cosx_x3(a, value):
    return a.cos - a * a * a - value


def septic(a, value):
    # Example: ./models.py 2 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def playground(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos - value
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr - value
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
    # return a.sqr.sqrt - value
    # return a.exp.ln - value
    # return (a.exp + (a ** 2 - 4).exp).sqrt - value
    return (a.exp + (a ** 2 - 4).exp).ln - value


analyze(septic)
