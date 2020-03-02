# /opt/Nim/Nim/bin/nim compile horner.nim
# ./horner

import taylor

let n = 4
let h = 3.0
let jet = [-19.0'f64, 7.0'f64, -4.0'f64, 6.0'f64]

echo(n)
echo(h)
echo(jet)
echo(len(jet))
echo(t_horner(jet, h))

