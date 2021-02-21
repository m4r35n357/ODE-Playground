#!/usr/bin/env python3
"""
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from json import loads
from sys import stdin, stderr, argv
from math import sqrt, acos, sin, cos, log10, fabs
from Symplectic import Dual, Symplectic

class BhSymp(object):
    def __init__(self, a, μ2, e, lz, cc, r0, θ0, xh):
        self.a = a
        self.μ2 = μ2
        self.E = e
        self.L = lz
        self.a2 = a**2
        self.a2μ2 = self.a2 * μ2
        self.aE = a * e
        self.aL = a * lz
        self.K = cc + (lz - self.aE)**2
        self.t = 0.0
        self.r = Dual.get(r0).var
        self.θ = Dual.get((90.0 - θ0) * acos(-1.0) / 180.0).var
        self.φ = 0.0
        self.cross = xh
        self.refresh()
        self.ur = - sqrt(self.R.val if self.R.val >= 0.0 else - self.R.val)
        self.uθ = - sqrt(self.Θ.val if self.Θ.val >= 0.0 else - self.Θ.val)

    def refresh(self):
        r2 = self.r**2
        self.ra2 = r2 + self.a2
        P = self.ra2 * self.E - self.aL
        self.Δ = self.ra2 - 2.0 * self.r
        self.R = P**2 - self.Δ * (self.μ2 * r2 + self.K)
        self.sin2θ = self.θ.sin**2
        cos2θ = 1.0 - self.sin2θ
        T = self.aE * self.sin2θ - self.L
        self.Θ = self.K - self.a2μ2 * cos2θ - T**2 / self.sin2θ
        P_Δ = P.val / self.Δ.val
        self.Σ = r2.val + self.a2 * cos2θ.val
        self.ut = P_Δ * self.ra2.val - T.val * self.a
        self.uφ = P_Δ * self.a - T.val / self.sin2θ.val

    def four_V(self):
        ut = self.ut / self.Σ;
        ur = self.ur / self.Σ;
        uθ = self.uθ / self.Σ;
        uφ = self.uφ / self.Σ;
        return (self.sin2θ.val / self.Σ * (self.a * ut - self.ra2.val * uφ)**2 + self.Σ / self.Δ.val * ur**2
                + self.Σ * uθ**2 - self.Δ.val / self.Σ * (ut - self.a * self.sin2θ.val * uφ)**2)

    def q_update(self, c):
        self.t += c * self.ut
        self.r.val += c * self.ur
        self.θ.val += c * self.uθ
        self.φ += c * self.uφ
        self.refresh()

    def p_update(self, d):
        self.ur += 0.5 * d * self.R.der
        self.uθ += 0.5 * d * self.Θ.der

    def solve(self, method, h, start, end, tr):
        mino = τ = 0.0
        step = 0
        while (step < end) and (self.cross or self.Δ.val > 0.0):
            if τ >= start and step % tr == 0:
                self.plot(mino, τ)
            method()
            step += 1
            mino = step * h
            τ += h * self.Σ
        self.plot(mino, τ)

    def log_error(self, error, floor):
        return 10.0 * log10(fabs(error) if fabs(error) >= floor else floor)

    def plot(self, mino, τ):
        v4_error = self.log_error(self.μ2 + self.four_V(), 1e-36)
        R_error = self.log_error((self.r.val * self.r.val - self.R.val) / self.Σ**2, 1e-36)
        THETA_error = self.log_error((self.θ.val * self.θ.val - self.Θ.val) / self.Σ**2, 1e-36)
        ra = sqrt(self.ra2.val)
        print(f'{ra*self.θ.sin.val*cos(self.φ):.6e} {ra*self.θ.sin.val*sin(self.φ):.6e} {self.r.val*cos(self.θ.val):.6e} {mino:.6e} {v4_error:.6e} {R_error:.6e} {THETA_error:.6e}')


if __name__ == "__main__":
    #  Example: ./Bh3d.py initial-conditions.json  | ./filegraphics-pi.py initial-conditions.json
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = open(argv[1]).read() if len(argv) == 2 else stdin.read()
    ic = loads(input_data)['IC']
    print(input_data, file=stderr)
    bh = BhSymp(ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['cross'])
    step = ic['step']
    bh.solve(Symplectic(bh, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
