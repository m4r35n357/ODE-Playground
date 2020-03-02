#  No URL yet

# /opt/Nim/Nim/bin/nim compile -d:release --gc:none wb.nim
# ./wb 16 8 0.1 10000 1 0 0 2.0 >/tmp/data; echo "splot '/tmp/data' with lines" | gnuplot -p

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

var x = newSeq[float](n + 1)
var y = newSeq[float](n + 1)
var z = newSeq[float](n + 1)

var tx = newSeq[float](n)
var sx = newSeq[float](n)
var wa = newSeq[float](n)
wa[0] = a

t_output(x0, y0, z0, 0.0)
for step in 1..steps+1:
    x[0] = x0
    y[0] = y0
    z[0] = z0
    for k in 0..n-1:
        tx[k] = t_tan_sec2(tx, sx, x, k, hyp)[0]
        x[k + 1] = (y[k] - x[k]) / float(k + 1)
        y[k + 1] = - t_prod(z, tx, k) / float(k + 1)
        z[k + 1] = (- wa[k] + t_prod(x, y, k) + t_abs(y, k)) / float(k + 1)
    x0 = t_horner(x, h)
    y0 = t_horner(y, h)
    z0 = t_horner(z, h)
    t_output(x0, y0, z0, float64(step) * h)

