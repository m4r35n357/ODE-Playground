
from sys import stderr

def lorentz(a, value):
    # Example: ./models.py 0 .001 .999 1001 0 1e-12 1e-12 | ./plotMany.py 1 10 >/dev/null
    return 1 / (1 - a.sqr).sqrt - value

def cosx_x3(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a * a.sqr - value

def septic(a, value):
    # Example: ./models.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6) - value

def composite1(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a.sqr - 4.0).exp).ln - value

def composite2(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.sqr + (a.exp - 4).sqr).sqrt - value

def playground(a, value):
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
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
    return (a + 7) / (3.0 - a)

print(__name__ + " module loaded", file=stderr)
