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

var x0 = parseFloat(params[4])
var y0 = parseFloat(params[5])
var z0 = parseFloat(params[6])

let tmp = parseFloat(params[7])
let gamma = parseFloat(params[8])

var x = newSeq[float](n + 1)
var y = newSeq[float](n + 1)
var z = newSeq[float](n + 1)

var a = newSeq[float](n)
var b = newSeq[float](n)
var c = newSeq[float](n)

var alpha = newSeq[float](n)
alpha[0] = tmp

var w1 = newSeq[float](n)
w1[0] = 1.0

t_output(x0, y0, z0, 0.0)
for step in 1..steps+1:
    x[0] = x0
    y[0] = y0
    z[0] = z0
    for k in 0..n-1:
        let x2_1 = t_sqr(x, k) - w1[k]
        a[k] = z[k] + x2_1
        b[k] = 3.0 * z[k] - x2_1
        c[k] = alpha[k] + t_prod(x, y, k)
        x[k + 1] = (t_prod(y, a, k) + gamma * x[k]) / float(k + 1)
        y[k + 1] = (t_prod(x, b, k) + gamma * y[k]) / float(k + 1)
        z[k + 1] = - 2.0 * t_prod(z, c, k) / float(k + 1)
    x0 = t_horner(x, h)
    y0 = t_horner(y, h)
    z0 = t_horner(z, h)
    t_output(x0, y0, z0, float64(step) * h)

