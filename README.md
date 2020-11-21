## Quick Start

#### Debian/Ubuntu packages:
```
sudo apt install bc git build-essential musl-tools pkg-config mesa-utils-extra python3-tk python3-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11
```
#### Python 3 Packages (for plotting), please use a virtual environment!
IMPORTANT - If virtualenvwrapper is newly installed, you need a fresh login to set up its environment.
```
mkvirtualenv -p /usr/bin/python3 taylor
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
#### Run an ODE simulation (ODE call):

Runs a named simulation, and prints results to stdout

tsm-\*-\* (c executables) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | x, y, z output precision in decimal places
2 | order of Taylor Series
3 | time step size
4 | number of steps
5,6,7 | initial conditions, x0, y0, z0
8+ | ODE parameters

##### Run & plot (3D plot using pi3d):
```
./tsm-thomas-dbg 6 10 0.1 30000 1 0 0 .185 | ./plot3d.py
```
##### Run & plot (animated matplotlib graph):
```
./tsm-thomas-dbg 6 10 0.1 30000 1 0 0 .185 | ./plotAnimated.py -5 5
```
##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-dbg 6 10 0.1 30000 1 0 0 .185 | gnuplot -p -e "set terminal wxt size 1200,900; splot '<cat' with lines"
```
##### Run & plot (2D gnuplot graph):
```
./tsm-thomas-dbg 6 10 0.1 30000 1 0 0 .185 >/tmp/$USER/data
gnuplot -p -e "set terminal wxt size 1200,900; plot '/tmp/$USER/data' using 4:1 with lines, '/tmp/$USER/data' using 4:2 with lines, '/tmp/$USER/data' using 4:3 with lines"
```

#### Bifurcation Diagrams:

Runs a simulation many times for different values of a single parameter, produces turning point data for plotting bifurcation diagrams in X, Y and Z, and saves plots to PNG files.

Optionally, skips initial transient by dropping first (datalines / value) results (10 is usually a good value).

bifurcation-scan (shell script) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | start of parameter range
2 | end of parameter range
3 | "transient skip" value (skip first lines / value, or 0)
4+ | ODE call with variable parameter replaced by ['$p']

##### Bifurcation Diagram (manual gnuplot graph):
```
./bifurcation-scan .1 .23 10 ./tsm-thomas-static 6 10 0.1 10000 1 0 0 '$p'
```
This produces three PNG files, one for each coordinate.
If you want to interact with actual plots (e.g. to read off parameter values for simulation), use a command like (for x):
```
gnuplot -p -e "set t wxt size 1350,800 background rgb 'grey85'; set grid back; plot '/tmp/bifurcationX' lt rgb 'dark-blue' w dots, '/tmp/bifurcationx' lt rgb 'dark-green' w dots"
```

#### Clean Numerical Simulation:

Runs a simulation twice, once with a "better" integrator, and shows the differences graphically.

cns (shell script) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | CNS function, Selects a better integrator for comparison, see below
2+ | ODE call

CNS function Parameter | Meaning
----------|-----------
step4 | The step size is reduced by a quarter
step8 | The step size is reduced by an eightth
order | The Taylor Series order is increased by two
both | The order is increased by one, and the step size by one half
both2 | The order is increased by two, and the step size by one quarter

##### CNS plot (matplotlib diff graph):
```
./cns both ./tsm-thomas-static 15 10 0.1 30000 1 0 0 .185
```

#### CNS Duration Scanning

Runs a simulation repeatedly with increasing order of integration, for each order showing the simulation time when the deviation threshold is exceeded.

cns-scan (shell script) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | Maximum order for Taylor integrator (minimum is 1)
2 | deviation threshold
3+ | ODE call

##### CNS duration vs. Simulation Order (gnuplot graph):
```
./cns-scan both 24 1 ./tsm-thomas-static 15 10 0.1 10000 1 0 0 .185 | gnuplot -p -e "plot '<cat' with lines"
```

#### Sensitivity to Initial Conditions:

Runs a simulation together with six additional ones (+- deviations in X, Y and Z axes)

ic (shell script) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | Initial separation between "original" trajectory and the extra ones
2+ | ODE call

##### Sensitivity to Initial Conditions (3D plot using pi3d):
```
./ic .001 ./tsm-thomas-static 6 10 0.1 30000 1 0 0 .185
```

For more background see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md
