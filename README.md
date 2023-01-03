
## STATUS: feature-complete
No planned new features, and no known bugs ;)

## OpenGL display

There is now an OpenGL executable built for each TSM model, and also separate n-body and black hole OpenGL programs.
Once you have built the project, use one of these yad UI scripts, or read the docs below to drive directly from the shell:
```
./ode-playground
./hamiltonian-playground
./blackhole-playground
./blackhole-generator
```
The n-body and black hole programs will most probably not be documented here, so in these cases the yad UI currently _is_ the documentation.

## yad UI
Form or list-based dialogue boxes for launching most common operations using yad

## Pure c99 (plus OpenGL)

Choice of Clang or GCC.
All console programs are written to and depend _only_ on the c99 standard and library.
External dependencies (OpenGL, FreeGLUT & GLEW) are only needed for plotters (*-gl).

All c floating point operations are executed in _long double_ precision.
This gives a choice of precision and performance on different platforms.

Platform | Long Double Implementation
----------|-----------
ix86 | 80 bit hardware float
x86-64 | 80 bit hardware float
armhf | 64 bit hardware float
aarch64 | 128 bit _software_ float

All scripts are in _Posix_ shell syntax.

## Quick Start

### Requirements - Debian/Ubuntu packages
```
sudo apt install bc git build-essential musl-tools gnuplot-x11 lcov freeglut3-dev glew-utils libglew-dev
```
Optional:
```
sudo apt install yad ffmpeg google-perftools libgoogle-perftools-dev
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground/pure_c
cd ODE-Playground
```

#### c Build (Clang by default)

```
make clean
make OR make CCC=gcc OR make CCC=clang OR make CCC=dbg
make test
```
There should be NO errors or warnings.

[UPDATE] kerr-image.c (not built automatically) shows warnings on arm64; it is 3rd party code

### ODE analysis with Taylor Integrators

These models use the arbitrary-order Taylor Series Method (TSM)

* Good selection of clients (models) included
* Investigate the validity of chaotic solutions against integrator order and step size
* Plot bifurcation diagrams, to find "interesting" parameter values to study

```
./ode-playground
```

### Hamiltonian analysis with Symplectic Integrators

2nd to 10th order Suzuki integrators, with a model to help in visualization
All models except N-Body use Dual Numbers for Automatic Differentiation

Examples:
* N-Body system
* Mass-spring system
* Central mass Newtonian system
* an analysis example to visualize the symplectic time-stepping sequence for each order of integration.

(The N-body example uses symplectic integration, but not dual numbers, because the differentiation is trivial in this case.)

No formal documentation yet, see the yad and c files for example usage.

```
./hamiltonian-playground
```

### Spinning Black Hole (Kerr) orbits using Symplectic Integrators

This example uses a "pseudo-Hamiltonian" approach with Dual Numbers to solve Carter's first-order equations of motion, separated in r and theta using Mino time.

No formal documentation yet, see the yad and c files for example usage.
```
./blackhole-playground
```
There are also programs to create model parameters and initial conditions for bounded particle and light orbits
```
./blackhole-generator
```

## Grubby details

#### Running c Tests

The tests enforce a "redundant web of densely interrelated functionality that cannot exist in the presence of coding errors"!

**libad-test** (c executable)

Parameter | Meaning
----------|-----------
1 | Order
2 | X value
3 | Error limit (absolute)
4 | (Optional) verbosity: 0: summary (default), 1: list, 2: detail

The final parameter can be set to 0 (or left absent) for a summary, 1 for individual tests, or 2 for full detail of Taylor Series.
Depending on the x value, some tests might be skipped owing to domain restrictions on some of the functions involved.

```
./libad-test 20 .5 1e-15

Horner
 23   23.000
153  153.000
201  201.000

Taylor Series Method: x'=1  y'=0  z'=-1
+1.000000000000e+00 +1.000000000000e+00 +1.000000000000e+00 0.000000e+00 _ _ _ 0.000
+1.105170918076e+00 +1.000000000000e+00 +9.048374180360e-01 1.000000e-01 _ _ _ 0.000
+1.221402758160e+00 +1.000000000000e+00 +8.187307530780e-01 2.000000e-01 _ _ _ 0.000
+1.349858807576e+00 +1.000000000000e+00 +7.408182206817e-01 3.000000e-01 _ _ _ 0.000
+1.491824697641e+00 +1.000000000000e+00 +6.703200460356e-01 4.000000e-01 _ _ _ 0.001
+1.648721270700e+00 +1.000000000000e+00 +6.065306597126e-01 5.000000e-01 _ _ _ 0.001
+1.822118800391e+00 +1.000000000000e+00 +5.488116360940e-01 6.000000e-01 _ _ _ 0.001
+2.013752707470e+00 +1.000000000000e+00 +4.965853037914e-01 7.000000e-01 _ _ _ 0.001
+2.225540928492e+00 +1.000000000000e+00 +4.493289641172e-01 8.000000e-01 _ _ _ 0.001
+2.459603111157e+00 +1.000000000000e+00 +4.065696597406e-01 9.000000e-01 _ _ _ 0.001
+2.718281828459e+00 +1.000000000000e+00 +3.678794411714e-01 1.000000e+00 _ _ _ 0.001
Check: e^1  e^0  e^-1
+2.718281828459e+00 +1.000000000000e+00 +3.678794411714e-01 1.000000e+00 _ _ _ 0.000

Recurrence Relations: x = 0.5
Total: 44, PASSED 44
```

**libdual-test** (c executable)

Parameter | Meaning
----------|-----------
1 | X value
2 | Error limit (absolute)
3 | (Optional) verbosity: 0: summary (default), 1: list, 2: detail

```
./libdual-test .5 1e-15

Dual Numbers: x = 0.5
Total: 40, PASSED 40

for i in .5 0 -.5; do ./libad-test 10 $i 1e-15; ./libdual-test $i 1e-15; done
```
#### Find examples for ODE parameters and many other things:
Useful commands are frequently added to the comments in source headings.
```
grep Example *
```
## Solving and Plotting ODEs
This use case only involves calling the "t-functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.

Where CPU timings are given, they are made on a Raspberry Pi 400, mildly overclocked to 2100MHz, and writing output to a tmpfs file.

#### Run a basic ODE simulation (ODE call):

Runs a named simulation, and prints results to stdout.
Each output line consists of a column each for x, y, z, t, followed by three turning point tags for generating bifurcation diagrams, and cumulative CPU usage.

**tsm-model-type-gl** (c executables)

Parameter | Meaning
----------|-----------
1 | Length of track(s)
2 | order of Taylor Series
3 | time step
4 | number of steps
5,6,7 | initial conditions, x0,y0,z0
8+ | Model parameters

##### Run & plot (OpenGL plot):
```
./tsm-thomas-gl 2000 10 0.1 30000 1 0 0 .185 >/tmp/$USER/data
```

**tsm-model-type** (c executables)

Parameter | Meaning
----------|-----------
1 | x,y,z output decimal places (0 for full precision binary hex)
2 | order of Taylor Series
3 | time step
4 | number of steps
5,6,7 | initial conditions, x0,y0,z0
8+ | Model parameters

##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-std 6 10 0.1 30000 1 0 0 .185 >/tmp/$USER/data
 gnuplot -p << EOF
set xyplane 0
set view 54.73561,135
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'
splot '/tmp/$USER/data' with lines
EOF
```
##### Run & plot (2D gnuplot graph):
```
./tsm-lorenz-std 6 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/$USER/data

 gnuplot -p << EOF
set terminal wxt size 1200,900
plot '/tmp/$USER/data' using 4:1 with lines, '' u 4:2 w l, '' u 4:3 w l
EOF
```
It should be possible to send output directly to gnuplot via a pipe, but many versions segfault when reading stdin so I now specify a temporary file instead.

### Bifurcation Diagrams:

This script runs a simulation many times for different values of a single parameter, and uses turning point tags in the ODE simulation output for plotting bifurcation diagrams in X, Y and Z, and saves plots to PNG files.
Because the simulations carry second derivatives in the coordinate jets, we can plot maxima and minima in different colours!
Optionally, we can skip the initial transient by dropping the first (datalines / value) results (10 is usually a good value).

**bifurcation-scan** (shell script)

Parameter | Meaning
----------|-----------
1 | start of parameter range
2 | end of parameter range
3 | "transient skip" value; skip first (lines / value), or 0
4+ | ODE call with variable parameter replaced by ['$p']

The general idea is to replace one of the model parameters with the string '$p' (including quotes!).

#### Bifurcation Diagram (manual gnuplot graph):

A fourth-order integrator is sufficient for bifurcation diagrams and will run faster; for this scenario we only care about transitions into and out of chaos, not accuracy within the chaotic regions.
Progress output is sent to /dev/null in the examples below for brevity, but is useful in most situations.
```
time -p ./bifurcation-scan .1 .225 10 ./tsm-thomas-std 6 4 0.1 10000 1 0 0 '$p' >/dev/null
Bifurcation Diagrams: [.1 .225 10 ./tsm-thomas-std 6 4 0.1 10000 1 0 0 $p]
real 194.46
user 195.12
sys 116.62
```
This produces three PNG files, one for each coordinate.
You can see them using any image viewer e.g. ImageMagick:
```
display /tmp/$USER/X.png
```
If you want to interact with actual plots (e.g. to read off parameter values for simulation), use a command like (for the x coordinate):
```
 gnuplot -p << EOF
set t wxt size 1350,800 background rgb 'grey85'
set grid back
plot '/tmp/$USER/bifurcationX' lt rgb 'dark-blue' with dots, '/tmp/$USER/bifurcationx' lt rgb 'dark-green' w d
EOF
```

#### Clean Numerical Simulation:

In a chaotic system, accuracy can only be maintained for a finite simulation time.
This script runs a given simulation twice, the second time with a "better" integrator, and shows the differences graphically.
The maximum achievable clean simulation time is determined by machine precision alone.
Step size and order can affect the efficiency of the calculations, but not the ultimate accuracy!

**cns** (shell script)

Parameter | Meaning
----------|-----------
1 | CNS function, Selects a better integrator for comparison, see below
2 | deviation threshold
3+ | ODE call

CNS function | Meaning
----------|-----------
step2 | The step size is halved (this is  now the _only_ "better" integrator!)
nosim | User-defined comparison between /tmp/$USER/dataA and /tmp/$USER/dataB

NOTE: the CPU time required for a clean simulation is of order:
```
O(required precision * order^2 * clean simulation time)
```
Since (empirically, for these Lorenz parameters and timestep)
* clean simulation time is approximately proportional to the required precision
* clean simulation time is approximately proportional to the Taylor series order

the CPU requirement can also be seen as:
```
O(required precision^4) or O(order^4) or O(clean simulation time^4)
```
##### Make a CNS plot:

These results for higher orders (16+) require 128-bit precision, i.e. aarch64 long double (software).
In hardware 80-bit (x86-64) or 64-bit (armhf) floating point, the maximum clean simulation time will be correspondingly lower.
```
./cns step2 1.0 ./tsm-lorenz-std 6 4 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz-std 6 8 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz-std 6 12 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz-std 6 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz-std 6 20 .01 10000 -15.8 -17.48 35.64 10 28 8 3

./cns step2 1.0 ./tsm-lorenz-std 6 28 .01 10000 -15.8 -17.48 35.64 10 28 8 3
```
If you need to re-plot after closing gnuplot, either use the "nosim" argument, or:
```
 gnuplot -p << EOF
set key horizontal left
plot '/tmp/$USER/dataA' using 4:1 t 'x' with lines lc black, '' u 4:2 t 'y' w l lc black, '' u 4:3 t 'z' w l lc black, '/tmp/$USER/dataB' using 4:1 t 'x' with lines lc 'forest-green', '' u 4:2 t 'y' w l lc 'dark-yellow', '' u 4:3 t 'z' w l lc 'dark-turquoise'
EOF
```

#### CNS Duration Scanning (TSM only)

Runs a simulation repeatedly with increasing order of integration, for each order showing the simulation time when the deviation threshold is exceeded.
You can run this to determine the maximum _useful_ integrator order to use, for a given step size.

**cns-scan** (shell script)

Parameter | Meaning
----------|-----------
1 | Maximum order for Taylor integrator (minimum is 2)
2 | deviation threshold
3+ | ODE call (with "order" argument set to "_")

#### CNS duration vs. Simulation Order (gnuplot graph) for the given step size:

The following commands perform a scan, and plot the simulation time and cpu time as histograms against integrator order:
```
./cns-scan 32 1 ./tsm-lorenz-std 6 _ .01 10000 -15.8 -17.48 35.64 10 28 8 3  | tee /tmp/$USER/data

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

##Finally
For more background on the Taylor Series Method for solving ODEs, see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md
