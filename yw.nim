#  https://www.deviantart.com/chaoticatmospheres/art/Strange-Attractors-The-Yu-Wang-Attractor-376644146

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none yw.nim
# ./yw 16 10 .001 50000 1 0 0 10 40 2 2.5 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

import os
import strutils
import taylor

let params = commandLineParams()

let n = parseInt(params[1])
let h = parseFloat(params[2])
let steps = parseInt(params[3])

var x0 = parseFloat(params[4])
var y0 = parseFloat(params[5])
var z0 = parseFloat(params[6])

let a = parseFloat(params[7])
let b = parseFloat(params[8])
let c = parseFloat(params[9])
let d = parseFloat(params[10])

var x = t_jet_c(n + 1, x0)
var y = t_jet_c(n + 1, y0)
var z = t_jet_c(n + 1, z0)

var xy = newSeq[float](n)
var e_xy = newSeq[float](n)

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..n-1:
        xy[k] = t_prod(x, y, k)
        x[k + 1] = a * (y[k] - x[k]) / float(k + 1)
        y[k + 1] = (b * x[k] - c * t_prod(x, z, k)) / float(k + 1)
        z[k + 1] = (t_exp(e_xy, xy, k) - d * z[k]) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)

