
from sys import stderr

def lorentz(v):
    # Example: ./models.py BI .001 .999 1001 7 1e-12 1e-12 | ./plotMany.py 1 10 >/dev/null
    # Example: ./series_test.py 7 .001 .999 1001 | ./plotMany.py 1 10 >/dev/null
    return (1.0 - v**2)**-0.5

def schwartzschild(r, e=0.962250, p_r=0.0, l_z=4.0):  # no t or phi!
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (e**2 / (1 - 2 / r) - p_r**2 * (1 - 2 / r) - (l_z / r)**2) / 2

def cosx_x3(a):
    # Example: ./models.py BI -8 8 1001 7 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a**3

def septic(a):
    # Example: ./models.py BI -8 8 1001 7 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6)

def composite1(a):
    # Example: ./models.py BI -8 8 1001 7 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a**2 - 4).exp).ln

def composite2(a):
    # Example: ./models.py BI -8 8 1001 7 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a**2 + (a.exp - 4)**2)**0.5 - 4

def playground(a):
    # Example: ./models.py BI -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./series_test.py 7 -8 8 1001 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos
    # return (2 * a).cos - a.cos**2 + a.sin**2
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin**3
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos**3
    # return (a**2 + 1.0e-15)**0.5  # problem at a == 0.0
    # return (a**2 + 1.0e-6)**-0.5   # problem at a == 0.0
    # return a.exp.ln - a
    # return a.tan - a.sin / a.cos  #!
    # return a.tanh - a.sinh / a.cosh
    # return (a + 7) / (3.0 - a)  # problem at a == 3.0
    # return abs(a)**2
    # return (abs(a) + 1.0e-6)**-4
    # return a / a.exp
    return (a - 1)**2 / (a.cosh + 1).ln - 1

print(__name__ + " module loaded", file=stderr)
