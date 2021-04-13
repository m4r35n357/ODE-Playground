## NEWS: Symplectic integrators using dual numbers for automatic differentiation

Two examples; Newton orbits and Kerr (black hole) orbits.
See c files for example usage.  They can also be piped to ./plot3d.py

## Quick Start

#### Debian/Ubuntu packages:
```
sudo apt install bc git build-essential musl-tools pkg-config mesa-utils-extra python3-tk python3-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11
```
#### Python 3 Packages (for plotting), please use a virtual environment!
```
mkvirtualenv -p /usr/bin/python3 taylor
pip install matplotlib pillow pi3d
```
#### Build (MUSL by default):
`./build && echo OK` or `./build gcc && echo OK` or `./build clang && echo OK`

There should be NO errors or warnings.  [UPDATE: kerr-image.c shows warnings on arm64; it is 3rd party code]

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
./tsm-lorenz-dbg 6 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-dbg 6 10 0.1 30000 1 0 0 .185 | gnuplot -p -e "set terminal wxt size 1200,900; splot '<cat' with lines"
```
##### Run & plot (2D gnuplot graph):
```
./tsm-lorenz-dbg 6 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/$USER/data
gnuplot -p -e "set terminal wxt size 1200,900; plot '/tmp/$USER/data' using 4:1 with lines, '/tmp/$USER/data' using 4:2 with lines, '/tmp/$USER/data' using 4:3 with lines"
```

#### Bifurcation Diagrams:

Run  a simulation many times for different values of a single parameter, uses turning point markers in the ODE simulation output for plotting bifurcation diagrams in X, Y and Z, and saves plots to PNG files.
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
You can see them using any image viewer e.g. ImageMagick:
```
display /tmp/$USER/X.png
```
If you want to interact with actual plots (e.g. to read off parameter values for simulation), use a command like (for x):
```
gnuplot -p -e "set t wxt size 1350,800 background rgb 'grey85'; set grid back; plot '/tmp/$USER/bifurcationX' lt rgb 'dark-blue' w dots, '/tmp/$USER/bifurcationx' lt rgb 'dark-green' w dots"
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
step4 | The step size is quartered
step8 | The step size is eightthed
order | The order is increased by two
both | The order is increased by one, and the step size halved
both2 | The order is increased by two, and the step size quartered

##### CNS plot (matplotlib diff graph):
```
./cns both ./tsm-lorenz-static 6 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
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
./cns-scan both 28 1 ./tsm-lorenz-static 6 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3  | gnuplot -p -e "plot '<cat' with boxes"
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

## Interactive Function Analysis with Python
This use case involves calling the Series methods and operators from ad.py, and the functions from plotters.py, from a Python interpreter.
This is how I like to set things up:
```
$ echo '
from math import *
from ad import *
from plotters import *
' > ipython.py
$ ipython -i ipython.py
```

##### Higher level function analysis
Plotting function and derivatives, together with root and turning point analysis:
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: f = lambda a: (a.exp + (a.sqr - 4.0).exp).ln

In [2]: scan_s(f)
NT  x: -1.962e+00  δx: +4.123e-16  f: +1.332e-15  \ ROOT___ 4
NT  x: -1.312e+00  δx: -0.000e+00  f: +0.000e+00  / MIN_MAX 4
NT  x: -1.849e-02  δx: -2.714e-13  f: +2.662e-13  / ROOT___ 3

In [3]: mplot_s(f)
```
(matplotlib plot not shown!)

##### Square root of two
Here is a quick example of function inversion.
There is a choice of analysis (root finding) method:
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: newton_d(lambda x: x * x - 2.0, x0=1.0)
Out[1]: Result(method='NT', x=1.414213562373095, f=4.440892098500626e-16, δx=-1.570092458683775e-16, count=6, sense='_', mode='ROOT___')

In [2]: timeit(newton_d(lambda x: x * x - 2.0, x0=1.0))
25.4 µs ± 1.35 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)

In [3]: bisect_d(lambda x: x * x, xa=1.0, xb=2.0)
Out[3]: Result(method='BI', x=2.0, f=4.0, δx=0.0, count=101, sense='_', mode='ROOT___')

In [4]: timeit(bisect_d(lambda x: x * x, xa=1.0, xb=2.0))
342 µs ± 7.4 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```

##### Automatic Differentiation
Here we calculate _all_ the derivatives of a simple cubic in x, followed by its sensitivities to each parameter a, b, c.
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: a = Series.get(5, 3.0)

In [2]: b = Series.get(5, 5.0)

In [3]: c = Series.get(5, 7.0)

In [4]: x = Series.get(5, 2.0)

In [5]: print(a * x**3 - b * x**2 + c * x - 5)
+1.300e+01 +0.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [6]: print(a * x.var**3 - b * x.var**2 + c * x.var - 5)
+1.300e+01 +2.300e+01 +1.300e+01 +3.000e+00 +0.000e+00

In [7]: print(a.var * x**3 - b * x**2 + c * x - 5)
+1.300e+01 +8.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [8]: print(a * x**3 - b.var * x**2 + c * x - 5)
+1.300e+01 -4.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [9]: print(a * x**3 - b * x**2 + c.var * x - 5)
+1.300e+01 +2.000e+00 +0.000e+00 +0.000e+00 +0.000e+00
```

##Finally
For more background on the Taylor Series Method for solving ODEs, see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md
and the taylor-ode.h in that branch.
