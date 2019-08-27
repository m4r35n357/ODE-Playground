#!/usr/bin/env python3
# Example: ./double.py 8 0.01 10001 1 1 1 1 1 0 1 0 | ./plotPi2d.py
from sys import argv
from math import sin, cos, log10
from ad import Dual

# def h(self, qr, pr, pφ):  # NOTE: qφ is absent from Hamiltonian -> conserved quantity pφ
#     return (pr**2 + pφ**2 / qr**2) / (2.0 * self.m) - self.gm / qr

def h(g, l1, m1, l2, m2, qθ1, pθ1, qθ2, pθ2):  # the Hamiltonian
    return (l2**2 * m2 * pθ1**2 + l1**2 * (m1 + m2) * pθ2**2 - 2.0 * m2 * l1 * l2 * pθ1 * pθ2 * (qθ1 - qθ2).cos) \
           / (2.0 * l1**2 * l2**2 * m2 * (m1 + m2 * (qθ1 - qθ2).sin**2)) - (m1 + m2) * g * l1 * qθ1.cos - m2 * g * l2 * qθ2.cos

def main():
    g, n, δt, steps = 1.0, int(argv[1]), float(argv[2]), int(argv[3])
    l1, m1, l2, m2 = float(argv[4]), float(argv[5]), float(argv[6]), float(argv[7])
    qθ1, pθ1, qθ2, pθ2 = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])  # initial values
    qθ1_s, pθ1_s, qθ2_s, pθ2_s = [0.0] * (n + 1), [0.0] * (n + 1), [0.0] * (n + 1), [0.0] * (n + 1)  # the derivative "jets"
    qθ1_d, pθ1_d, qθ2_d, pθ2_d = Dual.get(qθ1), Dual.get(pθ1), Dual.get(qθ2), Dual.get(pθ2)  # create the dual numbers
    h0 = h(g, l1, m1, l2, m2, qθ1_d, pθ1_d, qθ2_d, pθ2_d).val
    for step in range(1, steps + 1):
        qθ1_s[0], pθ1_s[0], qθ2_s[0], pθ2_s[0] = qθ1, pθ1, qθ2, pθ2  #  Taylor Series method for the Hamiltonian EOMs
        for k in range(n):
            qθ1_s[k + 1] = h(g, l1, m1, l2, m2, qθ1_d, pθ1_d.var, qθ2_d, pθ2_d).der / (k + 1)
            pθ1_s[k + 1] = - h(g, l1, m1, l2, m2, qθ1_d.var, pθ1_d, qθ2_d, pθ2_d).der / (k + 1)
            qθ2_s[k + 1] = h(g, l1, m1, l2, m2, qθ1_d, pθ1_d, qθ2_d, pθ2_d.var).der / (k + 1)
            pθ2_s[k + 1] = - h(g, l1, m1, l2, m2, qθ1_d, pθ1_d, qθ2_d.var, pθ2_d).der / (k + 1)
        qθ1, pθ1, qθ2, pθ2 = qθ1_s[n], pθ1_s[n], qθ2_s[n], pθ2_s[n]  #  Horner's method
        for i in range(n - 1, -1, -1):
            qθ1 = qθ1 * δt + qθ1_s[i]
            pθ1 = pθ1 * δt + pθ1_s[i]
            qθ2 = qθ2 * δt + qθ2_s[i]
            pθ2 = pθ2 * δt + pθ2_s[i]
        qθ1_d, pθ1_d, qθ2_d, pθ2_d = Dual.get(qθ1), Dual.get(pθ1), Dual.get(qθ2), Dual.get(pθ2)  # recreate the dual numbers
        x1 = l1 * sin(qθ1)  #  convert angles to X-Y coordinates
        y1 = - l1 * cos(qθ1)
        x2 = x1 + l2 * sin(qθ2)
        y2 = y1 - l2 * cos(qθ2)
        error = abs(h(g, l1, m1, l2, m2, qθ1_d, pθ1_d, qθ2_d, pθ2_d).val - h0)
        print(f"{x1:.9e} {y1:.9e} {x2:.9e} {y2:.9e} {step*δt:.5e} {10.0*log10(error if error > 1.0e-18 else 1.0e-18):.9e}")

main()
