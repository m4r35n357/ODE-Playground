#!/usr/bin/env python3

from sys import argv
from gmpy2 import get_context, sin, cos, log10
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from ad import Dual
from taylor import to_mpfr


def hamiltonian(g, l1, m1, l2, m2, th1, pth1, th2, pth2):
    return (l2**2 * m2 * pth1.sqr + l1**2 * (m1 + m2) * pth2.sqr
            - 2 * m2 * l1 * l2 * pth1 * pth2 * (th1 - th2).cos) \
           / (2 * l1**2 * l2**2 * m2 * (m1 + m2 * (th1 - th2).sin.sqr)) \
           - (m1 + m2) * g * l1 * th1.cos - m2 * g * l2 * th2.cos


def main():
    n = int(argv[1])
    h = to_mpfr(argv[2])
    steps = int(argv[3])
    l1, m1, l2, m2 = to_mpfr(argv[4]), to_mpfr(argv[5]), to_mpfr(argv[6]), to_mpfr(argv[7])
    th1_0, pth1_0, th2_0, pth2_0 = to_mpfr(argv[8]), to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
    g = to_mpfr(1)
    #  the derivative "jets"
    th1 = [to_mpfr(0) for _ in range(n + 1)]
    pth1 = [to_mpfr(0) for _ in range(n + 1)]
    th2 = [to_mpfr(0) for _ in range(n + 1)]
    pth2 = [to_mpfr(0) for _ in range(n + 1)]
    #  the dual numbers
    th1_dual = Dual.get(th1_0)
    pth1_dual = Dual.get(pth1_0)
    th2_dual = Dual.get(th2_0)
    pth2_dual = Dual.get(pth2_0)
    h0 = hamiltonian(g, l1, m1, l2, m2, th1_dual, pth1_dual, th2_dual, pth2_dual).val
    for step in range(1, steps + 1):
        th1[0], pth1[0], th2[0], pth2[0] = th1_0, pth1_0, th2_0, pth2_0
        #  generate derivatives for the Hamiltonian EOMs
        dh_dth1 = hamiltonian(g, l1, m1, l2, m2, th1_dual.var, pth1_dual, th2_dual, pth2_dual).der
        dh_dpth1 = hamiltonian(g, l1, m1, l2, m2, th1_dual, pth1_dual.var, th2_dual, pth2_dual).der
        dh_dth2 = hamiltonian(g, l1, m1, l2, m2, th1_dual, pth1_dual, th2_dual.var, pth2_dual).der
        dh_dpth2 = hamiltonian(g, l1, m1, l2, m2, th1_dual, pth1_dual, th2_dual, pth2_dual.var).der
        #  Taylor Series method
        for k in range(n):
            th1[k + 1] = dh_dpth1 / (k + 1)
            pth1[k + 1] = - dh_dth1 / (k + 1)
            th2[k + 1] = dh_dpth2 / (k + 1)
            pth2[k + 1] = - dh_dth2 / (k + 1)
        #  Horner's method
        th1_0, pth1_0, th2_0, pth2_0 = th1[n], pth1[n], th2[n], pth2[n]
        for i in range(n - 1, -1, -1):
            th1_0 = th1_0 * h + th1[i]
            pth1_0 = pth1_0 * h + pth1[i]
            th2_0 = th2_0 * h + th2[i]
            pth2_0 = pth2_0 * h + pth2[i]
        #  the dual numbers
        th1_dual = Dual.get(th1_0)
        pth1_dual = Dual.get(pth1_0)
        th2_dual = Dual.get(th2_0)
        pth2_dual = Dual.get(pth2_0)
        #  convert angles to X-Y coordinates
        x1 = l1 * sin(th1_0)
        y1 = - l1 * cos(th1_0)
        x2 = x1 + l2 * sin(th2_0)
        y2 = y1 - l2 * cos(th2_0)
        error = abs(hamiltonian(g, l1, m1, l2, m2, th1_dual, pth1_dual, th2_dual, pth2_dual).val - h0)
        print(f"{x1:.9e} {y1:.9e} {x2:.9e} {y2:.9e} {step*h:.5e} {10*log10(error if error > 1.0e-18 else 1.0e-18):.9e}")


main()
