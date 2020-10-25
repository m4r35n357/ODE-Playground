## Quick Start

Build (MUSL by default):
```
./build
```
Run:
```
./tsm-lorenz-dbg 15 _ 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plot3d.py
./tsm-lorenz-dbg 15 _ 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
Bifurcation Diagram:
```
./bifurcation-scan 2.5 50 ./tsm-lorenz-static 15 _ 10 .1 10000 -15.8 -17.48 35.64 10 '$p' 8 3
gnuplot -p -e "set terminal wxt size 1350,800; set grid back; plot '<cat' with dots" </tmp/bifurcationX
```
Clean Numerical Simulation
```
./cns both ./tsm-thomas-static 15 _ 10 0.1 30000 1 0 0 .185
```
CNS duration vs. Simulation Order:
```
./cns-scan both 24 1 ./tsm-lorenz-static 15 _ 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | gnuplot -p -e "plot '<cat' with lines"
```
Sensitivity to Initial Conditions:
```
./ic .001 ./tsm-thomas-static 15 _ 10 0.1 30000 1 0 0 .185
```
