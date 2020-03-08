# https://en.wikipedia.org/wiki/R%C3%B6ssler_attractor

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none rossler.nim
# ./rossler 16 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

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
let tmp = parseFloat(params[8])
let c = parseFloat(params[9])

var x = t_jet_c(n + 1, x0)
var y = t_jet_c(n + 1, y0)
var z = t_jet_c(n + 1, z0)

var b = newSeq[float](n)
b[0] = tmp

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..n-1:
        x[k + 1] = - (y[k] + z[k]) / float(k + 1)
        y[k + 1] = (x[k] + a * y[k]) / float(k + 1)
        z[k + 1] = (b[k] + t_prod(x, z, k) - c * z[k]) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)

