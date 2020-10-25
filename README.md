## Quick Start

Build (MUSL by default):
```
./build
```
Run:
```
./tsm-thomas-dbg 15 _ 10 0.1 30000 1 0 0 .185 | ./plot3d.py
./tsm-thomas-dbg 15 _ 10 0.1 30000 1 0 0 .185 | ./plotAnimated.py -5 5
```
Bifurcation Diagram:
```
./bifurcation-scan .1 .25 ./tsm-thomas-static 15 _ 10 0.1 10000 1 0 0 '$p'
gnuplot -p -e "set terminal wxt size 1350,800; set grid back; plot '<cat' with dots" </tmp/bifurcationX
```
Clean Numerical Simulation
```
./cns both ./tsm-thomas-static 15 _ 10 0.1 30000 1 0 0 .185
```
CNS duration vs. Simulation Order:
```
./cns-scan both 24 1 ./tsm-thomas-static 15 _ 10 0.1 10000 1 0 0 '$p' | gnuplot -p -e "plot '<cat' with lines"
```
Sensitivity to Initial Conditions:
```
./ic .001 ./tsm-thomas-static 15 _ 10 0.1 30000 1 0 0 .185
```
For more background see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md