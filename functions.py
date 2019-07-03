
from sys import stderr

def lorentz(a):
    # Example: ./models.py 0 .001 .999 1001 13 1e-12 1e-12 | ./plotMany.py 1 10 >/dev/null
    return 1 / (1 - a.sqr).sqrt

def cosx_x3(a):
    # Example: ./models.py 1 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a * a.sqr

def septic(a):
    # Example: ./models.py 2 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6)

def composite1(a):
    # Example: ./models.py 2 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a.sqr - 4.0).exp).ln

def composite2(a):
    # Example: ./models.py 1 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.sqr + (a.exp - 4).sqr).sqrt

def playground(a):
    # Example: ./models.py 2 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr
    # return a.sqr.sqrt - abs(a)
    # return (a**2)**0.5 - abs(a)
    # return a.exp.ln - a
    # return a.tan - a.sin / a.cos
    # return a.tanh - a.sinh / a.cosh
    # return a.sin.asin
    # return a.cos.acos
    # return a.tan.atan
    # return (a + 7) / (3.0 - a)
    # return a**4
    return (-4)**a
    # return (2 + 3 * a)**(2 * a - 5)

print(__name__ + " module loaded", file=stderr)
