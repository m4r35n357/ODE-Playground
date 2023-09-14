
## STATUS: feature-complete

No planned new features, no known bugs, good test coverage, clean static analysis.

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

The programs can be built either with Clang or GCC.
All console programs are written to and depend _only_ on the c99 standard and library (strictly speaking, the WG14/N1256 _draft_ standard!).
External dependencies (OpenGL, FreeGLUT & GLEW) are only needed for the OpenGL plotters (*-gl).

All simulation floating point operations are executed in _long double_ precision.
This gives a choice of precision and performance on different platforms.

Platform | Long Double Implementation
----------|-----------
x86 / x86-64 | 80 bit hardware float
armhf | 64 bit hardware float
aarch64 | 128 bit _software_ float

These are converted to _float_ for use by OpenGL.

All scripts are in _Posix_ shell syntax.

## Quick Start

### Requirements - Debian/Ubuntu packages
```
sudo apt install bc git build-essential gnuplot-qt lcov freeglut3-dev glew-utils libglew-dev
```
Optional:
```
sudo apt install yad feh ffmpeg google-perftools libgoogle-perftools-dev clang-tools
```
### Download
```
git clone https://github.com/m4r35n357/ODE-Playground/pure_c
cd ODE-Playground
```

### c Build (Clang by default)

```
make clean && make CCC=gcc
make test
```
There should be NO errors or warnings.

### ODE analysis with Taylor Integrators

These models use the arbitrary-order Taylor Series Method (TSM)

* Good selection of clients (models) included
* Investigate the validity of chaotic solutions against integrator order and step size
* Plot bifurcation (chaos scanning) diagrams, to find "interesting" parameter values to study

There is a script that provides a dialogue-based user interface for most of the supplied models:
```
./ode-playground
```

### Hamiltonian analysis with Symplectic Integrators

2nd to 10th order Suzuki integrators, with a model to help in visualization.
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

### Running c Tests

The tests enforce a "redundant mesh of densely interrelated functionality that cannot exist in the presence of coding errors", or something like that!

**libad-test** (c executable)

Parameter | Meaning
----------|-----------
1 | Order
2 | X value
3 | Error limit (absolute)
4 | (Optional) verbosity: 0: summary (default), 1: list, 2: detail

The final parameter can be set to 0 (or left absent) for a summary, 1 to list individual tests, or 2 for full detail of Taylor Series coefficients.
Depending on the x value, some tests might be skipped owing to domain restrictions on some of the functions involved.

```
./libad-test 20 .5 1e-15

argc: 4, argv: [ ./libad-test 20 .5 1e-15 ]
Horner Summation ... OK
Taylor Series Method .......... OK
Recurrence Relations x = 0.5
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

argc: 3, argv: [ ./libdual-test .5 1e-15 ]
Dual Numbers x = 0.5
Total: 40, PASSED 40
```

### Code coverage

```
make clean && make CCC=cov
make coverage
```
see local HTML link at the end of the output.
Note that any of the programs that you run from this point onwards will add to the coverage output.

### Profiling examples

**profile** (shell script)

Parameter | Meaning
----------|-----------
1  | gcc|gpt
2+ | ODE call (-std or -gl)

#### Google performance tools (requires google-perftools & libgoogle-perftools-dev packages on Debian)
```
./profile gpt ./tsm-lorenz-std  6 16 .01 1000000  -15.8 -17.48 35.64  10 28 8 3
```

#### GCC gprof
```
./profile gcc ./tsm-lorenz-std  6 16 .01 1000000  -15.8 -17.48 35.64  10 28 8 3
```

### Static Analysis (requires clang-tools package on Debian)
```
make clean && scan-build make
...
scan-build: No bugs found.
```

### Development

There is an optional git pre-commit script that you can use automatically by copying:
```
cp pre-commit .git/hooks
```
or it can be run manually:
```
./pre-commit
```
It does a clean build, runs tests, and performs basic sanity checks on key executables.

## Running the programs

### Solving and Plotting ODEs
This use case only involves calling the "t_functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors and http://www.atomosyd.net/spip.php?rubrique5.

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

#### Run & plot (OpenGL plot):
```
./tsm-thomas-gl 5000 8 0.100 10000 1.0 0.0 0.0 0.185
argc: 9, argv: [ ./tsm-thomas-gl 5000 8 0.100 10000 1.0 0.0 0.0 0.185 ]

  OpenGL: 3.1 Mesa 19.3.2
FreeGLUT: 20801
    GLEW: 2.1.0
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

#### Run & plot (3D gnuplot graph):
```
./plot ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185
args: 10, [ ./plot ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185 ]
argc: 9, argv: [ ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185 ]
```

### Bifurcation (chaos scanning) Diagrams:

This script runs a simulation many times for different values of a single parameter, and uses turning point tags in the ODE simulation output for plotting bifurcation diagrams in X, Y and Z, and saves plots to PNG files.
Because the simulations carry second derivatives in the coordinate jets, we can plot maxima and minima in different colours!
Optionally, we can skip the initial transient by dropping the first (datalines / value) results (10 is usually a good value).

**chaos-scan** (shell script)

Parameter | Meaning
----------|-----------
1 | start of parameter range
2 | end of parameter range
3 | number of valid samples at end of raw data
4 | save snapshots ? [ no-snaps | snaps ]
5+ | ODE call with variable parameter replaced by ['$p']

The general idea is to replace one of the model parameters with the string '$p' (including quotes!).

#### Chaos Scan (manual gnuplot graph):

A fourth-order integrator is sufficient for bifurcation diagrams and will run faster; for this scenario we only care about transitions into and out of chaos, not accuracy within the chaotic regions.
Progress output (current parameter value) is produced while the script is running.
```
time -p ./chaos-scan 0.177 0.198 10 ./tsm-thomas-std 6 4 0.100 20000 1.0 0.0 0.0 '$p'
args: 13, [ ./chaos-scan 0.177 0.198 10 ./tsm-thomas-std 6 4 0.100 20000 1.0 0.0 0.0 $p ]
real 92.67                                                  
user 95.45
sys 55.47
```
Timings are from a Raspberry Pi 4 clocked to 1800 MHz.
As well as three "interactive plots" (mouse hover reads out plot coordinates), this produces three PNG files, one each for X, Y, Z.
You can see them later using any image viewer e.g. ImageMagick (without coordinate readouts of course):
```
display /tmp/$USER/X.png
```
If you choose to save snapshots, you can view them later in a slideshow program of your choice, for example feh:
```
feh --slideshow-delay .1 /tmp/$USER/snapshots
```
or make a video:
```
ffmpeg -framerate 5 -pattern_type glob -i "/tmp/$USER/snapshots/*.png" -c:v libx264 -pix_fmt yuv420p ode.mp4
```
To suppress autoscaling you can provide axis limits to the chaos-scan script via environment variables:
```
time -p X_START=-3 X_END=3 Y_START=-4 Y_END=4 Z_START=0 Z_END=2 ./chaos-scan 0.1 0.35 10000 snaps ./tsm-rf-std 6 4 0.010 50000 0.05 -0.05 0.3 '$p' .1
```
you only need to add the ones that you need, and the order doesn't matter.

### Clean Numerical Simulation:

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
#### Make a CNS plot

These results for higher orders (16+) require 128-bit precision, i.e. aarch64 long double (software).
In hardware 80-bit (x86-64) or 64-bit (armhf) floating point, the maximum clean simulation time will be correspondingly lower.
```
./cns step2 1.0 ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185
args: 12, [ ./cns step2 1.0 ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185 ]
argc: 9, argv: [ ./tsm-thomas-std 6 8 .050000000 20000 1.0 0.0 0.0 0.185 ]
argc: 9, argv: [ ./tsm-thomas-std 6 8 0.100 10000 1.0 0.0 0.0 0.185 ]
argc: 4, argv: [ ./divergence /tmp/pi/dataA /tmp/pi/dataB 1.0 ]
```
If you need to re-plot after closing gnuplot, use the "nosim" argument.

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

The following command performs a scan, and plots the clean simulation time and cpu time as histograms against integrator order from 2 to 32:
```
./cns-scan 32 1.0 ./tsm-thomas-std 6 _ 0.100 10000 1.0 0.0 0.0 0.185
args: 12, [ ./cns-scan 32 1.0 ./tsm-thomas-std 6 _ 0.100 10000 1.0 0.0 0.0 0.185 ]
<output not shown>
```

##Finally
For more background on the Taylor Series Method for solving ODEs, see the old README:
https://github.com/m4r35n357/ODE-Playground/blob/master/README.md
