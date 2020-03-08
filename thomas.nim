#  https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none thomas.nim
# ./thomas 9 8 0.1 30000 1 0 0 .19 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

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

let b = parseFloat(params[7])

var x = t_jet_c(n + 1, x0)
var y = t_jet_c(n + 1, y0)
var z = t_jet_c(n + 1, z0)

var sx = newSeq[float](n)
var cx = newSeq[float](n)
var sy = newSeq[float](n)
var cy = newSeq[float](n)
var sz = newSeq[float](n)
var cz = newSeq[float](n)

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..n-1:
        x[k + 1] = (t_sin_cos(sy, cy, y, k, trig)[0] - b * x[k]) / float(k + 1)
        y[k + 1] = (t_sin_cos(sz, cz, z, k, trig)[0] - b * y[k]) / float(k + 1)
        z[k + 1] = (t_sin_cos(sx, cx, x, k, trig)[0] - b * z[k]) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)

