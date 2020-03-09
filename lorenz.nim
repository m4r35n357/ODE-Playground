# https://en.wikipedia.org/wiki/Lorenz_system

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none lorenz.nim
# ./lorenz 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

import os
import strutils
import taylor

let params = commandLineParams()

let n = parseInt(params[1])
let h = parseFloat(params[2])
let steps = parseInt(params[3])

var x = t_jet_c(n + 1, parseFloat(params[4]))
var y = t_jet_c(n + 1, parseFloat(params[5]))
var z = t_jet_c(n + 1, parseFloat(params[6]))

let sigma = parseFloat(params[7])
let rho = parseFloat(params[8])
let beta = parseFloat(params[9]) / parseFloat(params[10])

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..<n:
        x[k + 1] = sigma * (y[k] - x[k]) / float(k + 1)
        y[k + 1] = (rho * x[k] - y[k] - t_prod(x, z, k)) / float(k + 1)
        z[k + 1] = (t_prod(x, y, k) - beta * z[k]) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)

