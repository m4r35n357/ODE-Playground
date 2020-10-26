## Quick Start

#### Debian/Ubuntu packages:
```
sudo apt install bc git build-essential musl-tools pkg-config mesa-utils-extra python3-tk python3-dev libmpc-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11
```
#### Python Packages (for plotting)
```
pip install matplotlib pillow pi3d
```
#### Build (MUSL by default):
```
./build
```
#### Find examples for ODE parameters and other things:
```
grep Example *
```
#### Run an ODE simulation:

tsm-\*-\* c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | x, y, z output precision in decimal places
2 | order of Taylor Series
3 | time step size
4 | number of steps
5,6,7 | initial conditions, x0, y0, z0
8+ | ODE parameters

#### Run & plot (3D plot using pi3d):
```
./tsm-thomas-dbg 15 10 0.1 30000 1 0 0 .185 | ./plot3d.py
```
#### Run & plot (animated matplotlib graph):
```
./tsm-thomas-dbg 15 10 0.1 30000 1 0 0 .185 | ./plotAnimated.py -5 5
```
#### Bifurcation Diagrams:

bifurcation-scan shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | start of parameter range
2 | end of parameter range
3 | "transient skip" value (skip first 1/value lines, or 0)
4+ | ODE parameters with variable parameter replaced by ['$p']

#### Bifurcation Diagram (gnuplot graph):
```
./bifurcation-scan .1 .25 10 ./tsm-thomas-static 15 10 0.1 10000 1 0 0 '$p'
gnuplot -p -e "set terminal wxt size 1350,800; set grid back; plot '<cat' with dots" </tmp/bifurcationX
```
#### Clean Numerical Simulation:

cns shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
CNS function | Selects a better integrator for comparison, see below

CNS function Parameter | Meaning
----------|-----------
step4 | The step size is reduced by a quarter
step8 | The step size is reduced by an eightth
order | The Taylor Series order is increased by two
both | The order is increased by one, and the step size by one half
both2 | The order is increased by two, and the step size by one quarter

#### CNS plot (matplotlib diff graph):
```
./cns both ./tsm-thomas-static 15 10 0.1 30000 1 0 0 .185
```
#### CNS Scan

cns-scan shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
CNS function | Selects a better integrator for comparison
deviation | threshold value

#### CNS duration vs. Simulation Order (gnuplot graph):
```
./cns-scan both 24 1 ./tsm-thomas-static 15 10 0.1 10000 1 0 0 '$p' | gnuplot -p -e "plot '<cat' with lines"
```
#### Sensitivity to Initial Conditions:

ic shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
separation | Initial separation between "original" trajectory and the additional ones

#### Sensitivity to Initial Conditions (3D plot using pi3d):
```
./ic .001 ./tsm-thomas-static 15 10 0.1 30000 1 0 0 .185
```

For more background see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md
