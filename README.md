
## NEWS: Basic dialogue box UI for many operations using yad

Once you have set up your virtual environment, you can access many features from the UI using:
```
./ode-playground
```
Built-in models:
* tsm-lorenz.py
* tsm-thomas.py
* tsm-halvorsen.py
* tsm-rf.py
* rk4-lorenz.py
* rk4-thomas.py
* rk4-halvorsen.py
* rk4-rf.py

## Quick Start

### Requirements - Debian/Ubuntu packages
```
sudo apt install bc git build-essential musl-tools pkg-config mesa-utils-extra python3-tk python3-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11 lcov gnuplot-x11
```
Optional:
```
sudo apt install yad ffmpeg
```

#### Python 3 Packages, please use a virtual environment!
There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python executables, which share output "formats").
In the example invocations given below, communication between the executable and plotting script uses a Unix pipe.
The dependencies are:
* matplotlib for 2D graphs and progressive graphs
* pi3d for 3D progressive trajectories
* gnuplot for static 3D trajectories

This is a typical Python virtual environment setup:
```
mkvirtualenv --python /usr/bin/python3 taylor
pip install matplotlib pillow pi3d pytest pytest-cov mutmut ipython
```

#### Running Python Tests
Testing Requirements
* pytest
* pytest_cov
* mutmut (optional)

Most of the code is covered several times over.
Tests have been mutation tested with mutmut.
```
pytest ad_test.py solver_test.py -v
```

#### Find examples for ODE parameters and many other things:
Useful commands are frequently added to the comments in source headings.
```
grep Example *
```
## Solving and Plotting ODEs
This use case only involves calling the "t_functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.

Where CPU timings are given, they are made on a Raspberry Pi 400, mildly overclocked to 2100MHz, and writing output to a tmpfs file.

#### Run a basic ODE simulation (ODE call):

Runs a named simulation, and prints results to stdout.
Each output line consists of a column each for x, y, z, t, followed by three turning point tags for generating bifurcation diagrams, and cumulative CPU usage.

**tsm-model.py** (Python script)

| Parameter | Meaning                      |
|-----------|------------------------------|
| 1         | x,y,z output decimal places  |
| 2         | order of Taylor Series       |
| 3         | time step                    |
| 4         | number of steps              |
| 5,6,7     | initial conditions, x0,y0,z0 |
| 8+        | Model parameters             |

**rk4-model.py** (Python script)

| Parameter | Meaning                                                                      |
|-----------|------------------------------------------------------------------------------|
| 1         | x,y,z output decimal places                                                  |
| 2         | only plot every ? lines of output (for smaller internal step than plot step) |
| 3         | time step                                                                    |
| 4         | number of steps                                                              |
| 5,6,7     | initial conditions, x0,y0,z0                                                 |
| 8+        | Model parameters                                                             |

##### Run & plot (3D plot using pi3d):
```
./tsm-thomas.py 6 10 0.1 30000 1 0 0 .185 | ./plot3d.py
```
##### Run & plot (animated matplotlib graph):
```
./tsm-lorenz.py 6 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
##### Run & plot (3D gnuplot graph):
```
./plot ./tsm-thomas.py 6 10 0.1 30000 1 0 0 .185
```
It should be possible to send output directly to gnuplot via a pipe, but many versions segfault when reading stdin so I now specify a temporary file instead.

#### Clean Numerical Simulation:

In a chaotic system, accuracy can only be maintained for a finite simulation time.
This script runs a given simulation twice, the second time with a "better" integrator, and shows the differences graphically.
The maximum achievable clean simulation time is determined by machine precision alone.
Step size and order can affect the efficiency of the calculations, but not the ultimate accuracy!

**cns** (shell script)

| Parameter | Meaning                                                             |
|-----------|---------------------------------------------------------------------|
| 1         | CNS function, Selects a better integrator for comparison, see below |
| 2         | deviation threshold                                                 |
| 3+        | ODE call                                                            |

| CNS function | Meaning                                                                |
|--------------|------------------------------------------------------------------------|
| step2        | The step size is halved (this is  now the _only_ "better" integrator!) |
| nosim        | User-defined comparison between /tmp/$USER/dataA and /tmp/$USER/dataB  |

##### Make a CNS plot (matplotlib diff graph):

Here are some comparisons between TSM and RK4 for roughly similar clean simulation times in each case.
Note that RK4 quickly becomes impractical because of excessive CPU usage, whereas TSM can stay clean up to even higher time values.
These specific results require 128-bit precision, i.e. aarch64 long double (software).
In hardware 80-bit (x86-64) or 64-bit (armhf) floating point, the maximum clean simulation time will be correspondingly lower.
```
./cns step2 1.0 ./rk4-lorenz.py 6 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
./cns step2 1.0 ./tsm-lorenz.py 6 4 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./rk4-lorenz.py 6 10 .001 100000 -15.8 -17.48 35.64 10 28 8 3
./cns step2 1.0 ./tsm-lorenz.py 6 8 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./rk4-lorenz.py 6 100 .0001 1000000 -15.8 -17.48 35.64 10 28 8 3
./cns step2 1.0 ./tsm-lorenz.py 6 12 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./rk4-lorenz.py 6 1000 .00001 10000000 -15.8 -17.48 35.64 10 28 8 3 
./cns step2 1.0 ./tsm-lorenz.py 6 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./rk4-lorenz.py 6 10000 .000001 100000000 -15.8 -17.48 35.64 10 28 8 3 
./cns step2 1.0 ./tsm-lorenz.py 6 20 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz.py 6 28 .01 10000 -15.8 -17.48 35.64 10 28 8 3
```
If you need to re-plot after closing gnuplot, either use the "nosim" argument, or:
```
 gnuplot -p << EOF                                                             
set key horizontal left
plot '/tmp/$USER/dataA' using 4:1 t 'xA' with lines lc black, '' u 4:2 t 'yA' w l lc black, '' u 4:3 t 'zA' w l lc black, '/tmp/$USER/dataB' using 4:1 t 'xB' with lines lc 'red', '' u 4:2 t 'yB' w l lc 'green', '' u 4:3 t 'zB' w l lc 'blue'
EOF
```

#### CNS Duration Scanning (TSM only)

Runs a simulation repeatedly with increasing order of integration, for each order showing the simulation time when the deviation threshold is exceeded.
You can run this to determine the maximum _useful_ integrator order to use, for a given step size.

**cns-scan** (shell script) 

| Parameter | Meaning                                            |
|-----------|----------------------------------------------------|
| 1         | Maximum order for Taylor integrator (minimum is 2) |
| 2         | deviation threshold                                |
| 3+        | ODE call (with "order" argument set to "_")        |

#### CNS duration vs. Simulation Order (gnuplot graph) for the given step size:

The following commands perform a scan, and plot the simulation time and cpu time as histograms against integrator order:
```
./cns-scan 32 1 ./tsm-lorenz.py 6 _ .01 10000 -15.8 -17.48 35.64 10 28 8 3  | tee /tmp/$USER/data

 gnuplot -p << EOF
set key left
set ytics nomirror
set y2tics
set xlabel 'Taylor Series Order'
set ylabel 'CNS Time, model units'
set y2label 'CPU Time, seconds'
plot '/tmp/$USER/data' using 1:2 axes x1y1 title 'CNS' with boxes, '' u 1:3 axes x1y2 t 'CPU' w boxes
EOF
```

#### Sensitivity to Initial Conditions (3D plot using pi3d):

Runs a simulation together with six additional ones (+- deviations in X, Y and Z axes) and plots directly to Pi3D.

**ic** (shell script)

| Parameter | Meaning                                                             |
|-----------|---------------------------------------------------------------------|
| 1         | Initial separation between "original" trajectory and the extra ones |
| 2         | Precision in decimal places ("scale" variable in bc)                |
| 3+        | ODE call                                                            |

The simulation is run seven times in parallel processes, the original along with each perturbed x, y, z.
```
./ic .001 32 ./tsm-thomas.py 6 10 0.1 30000 1 0 0 .185
```
```
./ic .001 32 ./tsm-lorenz.py 6 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 2>/dev/null
```
(3D plots not shown!)

### Taylor Series & Dual Numbers Calculator - Interactive Function Analysis with Python
This use case involves calling the Dual and Series methods and operators from ad.py, together with the functions from analysis.py, from a Python interpreter.
This is how I like to set things up (pip installed ipython recommended!):
```
$ echo '
from math import *
from ad import *
from analysis import *
' > ipython.py
$ ipython3 -i ipython.py
```

### Higher level function analysis (via the Python console or ipython)

Plotting function and derivatives, together with root and turning point analysis.
Scanning and plotting range is -8.0 to 8.0 by default.
Roots are found via zeros of f[0], turning points by zeros of f[1], and inflections via zeros of f[2].

* scan_d() finds roots only.
* mplot_d() plots the function and its first derivative.
* scan_s() finds roots, turning points and inflections.
* mplot_s() plots the function and its first 12 derivatives by default; 

Partly to verify my own implementation of the Taylor recurrence rules, I have added a demonstration of using series arithmetic to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Of course, it is more generally useful for finding inverse values (where real solutions exist) of complicated functions, not just their roots.

The "higher-level" ad_functions (or the Series class in Python) manipulate entire  Taylor series "jets" at once, so are only useful for univariate functions.
In Python there is an overloaded operator "~" which extracts the actual derivative values from the taylor series coefficents, as well as additional operators for (negation and **).
The \*\* (power) operator caters for f(x)^a, a^f(x) and f1(x)^f2(x), subject to domain limitations on f(x).
There are also functions for (matching the lower-level t_functions):
* abs
* sqr
* sqrt
* exp
* sin_cos (trig & hyp)
* tan_sec2 (trig & hyp)
* pwr (f(x)^a, where a is a scalar)
* ln
* asin, acos, atan (trig & hyp)

Using these "higher level" functions, Newton's method is implemented trivially, but I have also provided an implementation of the bisection method for comparison.

The analysis.py script enables the analysis of a "model" function's derivatives along with the value itself, across a range of the input variable.
That file contains a selection of example invocations in the comments.
Optionally it will analyse that range for roots, extrema and inflection points using the lower order derivatives.
Here is a seventh-degree polynomial model by way of example, and a trigonometric identity as a check of the sin and sqr functions:
```python
def septic(a):
    return (a + 7) * (5 + a) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6)

def trig(a):
    return (3 * a).sin - 3 * a.sin + 4 * a.sin**3
```
Operators on or between jets and numbers are implemented by overloading, and the functions are implemented as _properties_ of a jet.
These two approaches produce a fairly readable format for coding the models.
The value parameter allows simple calculation of function inverse values, as well as roots.

```
$ ipython3 -i ipython.py 
Python 3.9.2 (default, Feb 28 2021, 17:03:44) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.20.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
analysis module loaded

In [1]: f = lambda a: (a.exp + (a.sqr - 4.0).exp).ln

In [2]: scan_s(f)
NT  x: -1.962e+00  δx: +4.123e-16  f: +1.332e-15  \ ROOT___ 4
NT  x: -1.312e+00  δx: -0.000e+00  f: +0.000e+00  / MIN_MAX 4
NT  x: -1.849e-02  δx: -2.714e-13  f: +2.662e-13  / ROOT___ 3

In [3]: mplot_s(f)
```
(matplotlib plot not shown!)

### Square root of two
Here is a quick example of using root finding to invert a function.
The problem is set up as: [function] - [target value] = 0.
There is a choice of analysis (root finding) method:
You need to provide an initial estimate for the Newton solver, or an initial bracket for Bisection.
```
$ ipython3 -i ipython.py 
Python 3.9.2 (default, Feb 28 2021, 17:03:44) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.20.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
analysis module loaded

In [1]: newton_d(lambda x: x * x - 2.0, x0=1.0)
Out[1]: Result(method='NT', x=1.414213562373095, f=4.440892098500626e-16, δx=-1.570092458683775e-16, count=6, sense='_', mode='ROOT___')

In [2]: timeit(newton_d(lambda x: x * x - 2.0, x0=1.0))
19.2 µs ± 118 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)

In [3]: bisect_d(lambda x: x * x - 2.0, xa=1.0, xb=2.0)
Out[3]: Result(method='BI', x=1.414213562372879, f=-6.108447081487611e-13, δx=4.547473508864641e-13, count=41, sense='_', mode='ROOT___')

In [4]: timeit(bisect_d(lambda x: x * x - 2.0, xa=1.0, xb=2.0))
140 µs ± 691 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
```

### Automatic Differentiation
Here we calculate the value of a simple cubic in x(5), _all_ its derivatives wrt x(6), followed by the sensitivities to each parameter a(7), b(8), c(9).
Note that the function value is the same in each case, as you would expect!
```
$ ipython3 -i ipython.py 
Python 3.9.2 (default, Feb 28 2021, 17:03:44) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.20.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
analysis module loaded

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
