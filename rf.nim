#  https://en.wikipedia.org/wiki/Rabinovich%E2%80%93Fabrikant_equations

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none rf.nim
# ./rf 16 10 .01 50000 .05 -.05 .3 .28713 .1 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

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

let tmp = parseFloat(params[7])
let gamma = parseFloat(params[8])

var a = newSeq[float](n)
var b = newSeq[float](n)
var c = newSeq[float](n)

var alpha = t_jet_c(n, tmp)
var w1 = t_jet_c(n, 1.0)

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..<n:
        let x2_1 = t_sqr(x, k) - w1[k]
        a[k] = z[k] + x2_1
        b[k] = 3.0 * z[k] - x2_1
        c[k] = alpha[k] + t_prod(x, y, k)
        x[k + 1] = (t_prod(y, a, k) + gamma * x[k]) / float(k + 1)
        y[k + 1] = (t_prod(x, b, k) + gamma * y[k]) / float(k + 1)
        z[k + 1] = - 2.0 * t_prod(z, c, k) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)
