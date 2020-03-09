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

var x = t_jet_c(n + 1, parseFloat(params[4]))
var y = t_jet_c(n + 1, parseFloat(params[5]))
var z = t_jet_c(n + 1, parseFloat(params[6]))

let a = parseFloat(params[7])

var tx = newSeq[float](n)
var sx = newSeq[float](n)

var wa = t_jet_c(n, a)

t_output(x[0], y[0], z[0], 0.0)
for step in 1..steps+1:
    for k in 0..<n:
        tx[k] = t_tan_sec2(tx, sx, x, k, hyp)[0]
        x[k + 1] = (y[k] - x[k]) / float(k + 1)
        y[k + 1] = - t_prod(z, tx, k) / float(k + 1)
        z[k + 1] = (- wa[k] + t_prod(x, y, k) + t_abs(y, k)) / float(k + 1)
    x[0] = t_horner(x, h)
    y[0] = t_horner(y, h)
    z[0] = t_horner(z, h)
    t_output(x[0], y[0], z[0], float(step) * h)

