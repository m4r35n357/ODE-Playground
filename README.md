## Background

This project is mainly a collection of programs in c for evolving systems of ODEs using the Taylor Series Method (TSM), a rather old but poorly acknowledged technique based on forward mode Automatic Differentiation (AD).
TSM is a procedure for integrating ODEs using Taylor Series of arbitrary order, using recurrence relations between time derivatives of increasing order.
It is therefore (or should be, if it were better known!) a serious competitor to the fourth-order RK4 for the vast majority of cases.

The code itself is tiny (as are the dynamic executables) and has been partly developed on, and is suitable for running on, a Raspberry Pi 4 computer.

For the uninitiated, here is a review of the TSM itself and its history of repeated "discovery" and re-branding.
https://arxiv.org/abs/1111.7149

My work was inspired by parts of the following paper and its associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

Their software (called "Taylor") uses code generation to produce the recurrences and code to drive them, and it inspired me to try coding something more direct by hand.
These programs are the result.

Finally, more general resources on automatic differentiation can be found at the following portal: http://www.autodiff.org/ (see specifically the "Applications" and "Tools" sections).

## Taylor Series ODE Solvers (c/MPFR and Python)

My primary aim was to be able to solve coupled nonlinear equations and investigate chaotic systems, without relying on "black-box" ODE solvers.
The resulting c code takes the form of a small (<150 loc) Taylor Series "library", and the model-specific ODE simulators are tiny client programs to this, typically 15-25 loc each.
The header file taylor-ode.h contains a terse but complete description of the Taylor Series Method as implemented here, together with derivations of the recurrences that enable analysis of complex composed functions.

The recurrence rules (the "t-functions" in c and Python) are the key to calculating high order derivatives accurately, without needing finite differences.
They generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order values.
The functions provided cover the basic algebraic operations on Taylor Series (+ - * /), and also include several common functions:
* abs
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, where a is a scalar)
* ln

The recurrence relations used here are derived along the lines of (amongst other sources) http://www2.math.uni-wuppertal.de/wrswt/preprints/prep_05_4.pdf and http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).

There are also convenient factories for generating derivative "jets" of arbitrary order, and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, when properly called,  are all that is needed to solve systems of ODEs.
There is a fairly extensive collection of nonlinear ODE examples already implemented, in the file tsm.py, and in the tsm-\*-*.c files.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and many others.

## Accuracy of solutions
Rather than estimating local error using the usual Taylor Series method, I have provided a clean numerical simulation (CNS) shell script, that makes it easy to run a "better" simulation alongside.
The differences can be plotted for easy visual comparison.

## Scanning for chaos the simple way

There have been significant recent developments in automated testing for chaos, in particular from Gottwald & Melbourne http://www.maths.usyd.edu.au/u/gottwald/preprints/testforchaos_MPI.pdf and Wernecke, https://arxiv.org/abs/1605.05616.

These tend to be fairly involved techniques using various transforms and statistical methods against ensembles of nearby trajectories.
The method presented here is a very basic and literal approach based on deviations from a nominal trajectory, using a technique very similar to what is already used in the **Clean Numerical Simulation** shell script described below.

It was motivated by the first part of the Wernecke paper, but does not attempt the full method, nor use any statistical tools such as ensemble average, or correlation.
It should be seen as a solution fingerprinting tool for identifying regions of interest rather than a strict '0-1' test, which considering the fractal nature of the results (amongst other factors), is probably an unattainable objective.
As well as distinguishing between limit cycles and chaos, the method used here also identifies unbounded and converged solutions.

In summary, there are three main areas of application for the code:
* solving nonlinear ODEs
* scanning ODE parameters for detecting chaos in solutions
* checking global accuracy of solutions

## Build/Test Environment (Debian/Ubuntu/Raspbian)

OS level requirements:
```
sudo apt install bc git build-essential pkg-config mesa-utils-extra python3-tk python3-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground
cd ODE-Playground
```

#### Graph Plotting
There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python executables, which share output "formats").
In the example invocations given below, communication between the executable and plotting script uses a Unix pipe.
The dependencies are:
* matplotlib for 2D graphs and progressive graphs
* pi3d for 3D progressive trajectories
* gnuplot for static 3D trajectories

This is a typical Python virtual environment setup:
```
mkvirtualenv --python /usr/bin/python3 ad
pip install matplotlib pillow pi3d
```
Now you can just use it "in place" in your new virtual environment.

#### c Build (GCC or Clang)
```
$ ./build
$ ./build clang
```

##### Build command
```
./build && echo OK
```

## Solving and Plotting ODEs
This use case only involves calling the "t-functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm.py and tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.
To find some example invocations:
```
$ grep Example *.py tsm-*.c
```
Matplotlib progressive graph plotting:
```
$ ./tsm-lorenz-dbg NA NA 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```

tsm-\*-\* c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | Not used
2 | Not used
3 | order of Taylor Series (plot interval in RK4)
4 | step size
5 | number of steps
6,7, 8 | x0, y0, z0
9+ | ODE parameters

Since the RK4 Lorentz simulator (c only) is by definition fixed order, the "order" parameter is used to pass in a "plot interval" (e.g. 10 means plot only every 10th result).

#### 3D static trajectory plotting (gnuplot)

Write to a data file:
```
$ ./tsm-lorenz-dbg NA NA 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3  >/tmp/data
```
then plot it:
```
$ gnuplot -p -e "splot '/tmp/data' with lines"
```

#### 3D progressive trajectory plotting (pi3d)
```
$ ./tsm-lorenz-dbg NA NA 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3  | ./plot3d.py
```

## Scanning for chaos

This is now a bit easier to use, with a small collection of pre-supplied ODE models.

The chaos-scan shell script calls the ic shell script to generate six "nearby" trajectories over a range of a chosen ODE parameter.
These trajectories are then processed in the following order by chaos-distance.py to generate a colour-coded text summary, with embedded data suitable for plotting.
* Unbounded solutions are identified by length of the truncated results file
* Converged solutions are identified by the final separation being *smaller* than the initial separation
* Limit cycles are identified by *proportionality* of final to to initial separation over two different runs at different separations
* Chaotic solutions are identified by final separation being *independent* of the initial separation
* Remaining solutions are marked as UNCLASSIFIED for further investigation

chaos-scan shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
model | [ thomas lorenz rossler sprott halvorsen rf ]
start | Start of parameter range
end | End of parameter range
datalines | Number of time values to expect (for detecting unbounded solutions)
separation | Initial separation of the six main perturbed values
ratio | Relative separation of the six smaller deviation values

The script produces 100 lines of output, equally spaced between start and end.
The last three parameters have default values set (see script for details) and can usually be omitted:
```
$ time -p ./chaos-scan thomas .010 .250 2>/dev/null | tee /tmp/results
.010 UNCLASSIFIED value = 5.159e+01 3.624e+01 ratio = 1.423
.012400000000000      CHAOTIC value = 4.252e+01 4.482e+01 ratio = 0.949
.014800000000000      CHAOTIC value = 3.596e+01 3.539e+01 ratio = 1.016
.017200000000000 UNCLASSIFIED value = 3.604e+01 4.225e+01 ratio = 0.853
.019600000000000      CHAOTIC value = 3.545e+01 3.554e+01 ratio = 0.997
.022000000000000 UNCLASSIFIED value = 3.954e+01 3.481e+01 ratio = 1.136
.024400000000000      CHAOTIC value = 3.511e+01 3.388e+01 ratio = 1.036
.026800000000000      CHAOTIC value = 3.503e+01 3.688e+01 ratio = 0.950
.029200000000000 UNCLASSIFIED value = 3.278e+01 2.827e+01 ratio = 1.160
.031600000000000 UNCLASSIFIED value = 2.527e+01 2.814e+01 ratio = 0.898
.034000000000000      CHAOTIC value = 2.874e+01 2.850e+01 ratio = 1.008
.036400000000000      CHAOTIC value = 3.021e+01 2.856e+01 ratio = 1.058
.038800000000000      CHAOTIC value = 2.554e+01 2.623e+01 ratio = 0.974
.041200000000000 UNCLASSIFIED value = 2.315e+01 2.600e+01 ratio = 0.890
.043600000000000      CHAOTIC value = 2.526e+01 2.671e+01 ratio = 0.946
.046000000000000      CHAOTIC value = 2.429e+01 2.334e+01 ratio = 1.041
.048400000000000      CHAOTIC value = 2.309e+01 2.239e+01 ratio = 1.032
.050800000000000      CHAOTIC value = 2.226e+01 2.102e+01 ratio = 1.059
.053200000000000 UNCLASSIFIED value = 2.425e+01 2.133e+01 ratio = 1.137
.055600000000000      CHAOTIC value = 2.139e+01 2.197e+01 ratio = 0.974
.058000000000000      CHAOTIC value = 1.814e+01 1.971e+01 ratio = 0.921
.060400000000000      CHAOTIC value = 1.960e+01 1.980e+01 ratio = 0.990
.062800000000000      CHAOTIC value = 1.910e+01 1.851e+01 ratio = 1.032
.065200000000000      CHAOTIC value = 1.914e+01 1.875e+01 ratio = 1.021
.067600000000000      CHAOTIC value = 1.863e+01 1.911e+01 ratio = 0.975
.070000000000000      CHAOTIC value = 1.868e+01 1.849e+01 ratio = 1.010
.072400000000000 UNCLASSIFIED value = 1.662e+01 1.872e+01 ratio = 0.888
.074800000000000      CHAOTIC value = 1.862e+01 1.800e+01 ratio = 1.035
.077200000000000      CHAOTIC value = 1.735e+01 1.852e+01 ratio = 0.937
.079600000000000      CHAOTIC value = 1.801e+01 1.712e+01 ratio = 1.052
.082000000000000      CHAOTIC value = 1.846e+01 1.707e+01 ratio = 1.081
.084400000000000      CHAOTIC value = 1.820e+01 1.706e+01 ratio = 1.066
.086800000000000      CHAOTIC value = 1.710e+01 1.589e+01 ratio = 1.076
.089200000000000      CHAOTIC value = 1.583e+01 1.629e+01 ratio = 0.971
.091600000000000 UNCLASSIFIED value = 1.674e+01 1.387e+01 ratio = 1.207
.094000000000000 UNCLASSIFIED value = 1.480e+01 1.313e+01 ratio = 1.127
.096400000000000      CHAOTIC value = 1.244e+01 1.150e+01 ratio = 1.082
.098800000000000 UNCLASSIFIED value = 1.387e+01 1.254e+01 ratio = 1.106
.101200000000000      CHAOTIC value = 1.054e+01 1.055e+01 ratio = 0.998
.103600000000000 UNCLASSIFIED value = 7.774e+00 7.033e+00 ratio = 1.105
.106000000000000 UNCLASSIFIED value = 4.642e+00 5.167e+00 ratio = 0.898
.108400000000000  LIMIT-CYCLE value = 4.363e-03 4.403e-04 ratio = 9.909
.110800000000000  LIMIT-CYCLE value = 2.179e-04 2.178e-05 ratio = 10.002
.113200000000000  LIMIT-CYCLE value = 5.612e-04 5.612e-05 ratio = 10.001
.115600000000000  LIMIT-CYCLE value = 1.942e-04 1.942e-05 ratio = 10.000
.118000000000000  LIMIT-CYCLE value = 6.280e-05 6.280e-06 ratio = 10.000
.120400000000000  LIMIT-CYCLE value = 5.339e-05 5.339e-06 ratio = 10.000
.122800000000000 UNCLASSIFIED value = 9.322e+00 8.418e+00 ratio = 1.107
.125200000000000  LIMIT-CYCLE value = 1.006e-05 1.006e-06 ratio = 10.000
.127600000000000  LIMIT-CYCLE value = 5.592e-05 5.592e-06 ratio = 10.000
.130000000000000  LIMIT-CYCLE value = 1.088e-05 1.088e-06 ratio = 10.000
.132400000000000  LIMIT-CYCLE value = 1.134e-05 1.134e-06 ratio = 10.000
.134800000000000 UNCLASSIFIED value = 8.322e-01 5.878e-02 ratio = 14.159
.137200000000000  LIMIT-CYCLE value = 4.220e-03 4.209e-04 ratio = 10.027
.139600000000000  LIMIT-CYCLE value = 3.599e-05 3.599e-06 ratio = 10.000
.142000000000000  LIMIT-CYCLE value = 1.922e-05 1.922e-06 ratio = 10.000
.144400000000000  LIMIT-CYCLE value = 1.762e-05 1.762e-06 ratio = 10.000
.146800000000000      CHAOTIC value = 8.289e+00 8.878e+00 ratio = 0.934
.149200000000000 UNCLASSIFIED value = 2.429e-01 2.101e-02 ratio = 11.564
.151600000000000  LIMIT-CYCLE value = 5.227e-03 5.225e-04 ratio = 10.003
.154000000000000  LIMIT-CYCLE value = 1.341e-05 1.341e-06 ratio = 10.000
.156400000000000  LIMIT-CYCLE value = 3.124e-05 3.124e-06 ratio = 10.000
.158800000000000  LIMIT-CYCLE value = 1.194e-04 1.194e-05 ratio = 10.000
.161200000000000  LIMIT-CYCLE value = 4.972e-05 4.972e-06 ratio = 10.000
.163600000000000  LIMIT-CYCLE value = 6.120e-05 6.120e-06 ratio = 10.000
.166000000000000  LIMIT-CYCLE value = 1.125e-04 1.125e-05 ratio = 10.000
.168400000000000  LIMIT-CYCLE value = 2.569e-04 2.569e-05 ratio = 10.000
.170800000000000  LIMIT-CYCLE value = 1.014e-04 1.014e-05 ratio = 10.000
.173200000000000  LIMIT-CYCLE value = 6.226e-05 6.226e-06 ratio = 10.000
.175600000000000  LIMIT-CYCLE value = 4.336e-05 4.336e-06 ratio = 10.000
.178000000000000  LIMIT-CYCLE value = 3.234e-05 3.234e-06 ratio = 10.000
.180400000000000  LIMIT-CYCLE value = 2.605e-05 2.605e-06 ratio = 10.000
.182800000000000  LIMIT-CYCLE value = 1.971e-05 1.971e-06 ratio = 10.000
.185200000000000      CHAOTIC value = 9.171e+00 9.168e+00 ratio = 1.000
.187600000000000      CHAOTIC value = 9.152e+00 9.166e+00 ratio = 0.998
.190000000000000  LIMIT-CYCLE value = 1.935e-05 1.935e-06 ratio = 10.000
.192400000000000      CHAOTIC value = 9.135e+00 9.137e+00 ratio = 1.000
.194800000000000      CHAOTIC value = 9.107e+00 9.086e+00 ratio = 1.002
.197200000000000      CHAOTIC value = 9.025e+00 9.072e+00 ratio = 0.995
.199600000000000      CHAOTIC value = 6.379e+00 6.292e+00 ratio = 1.014
.202000000000000      CHAOTIC value = 6.198e+00 6.086e+00 ratio = 1.018
.204400000000000  LIMIT-CYCLE value = 2.556e-04 2.556e-05 ratio = 10.002
.206800000000000      CHAOTIC value = 5.964e+00 6.041e+00 ratio = 0.987
.209200000000000  LIMIT-CYCLE value = 2.751e-04 2.751e-05 ratio = 10.000
.211600000000000  LIMIT-CYCLE value = 1.283e-04 1.283e-05 ratio = 10.000
.214000000000000  LIMIT-CYCLE value = 4.019e-05 4.019e-06 ratio = 10.000
.216400000000000  LIMIT-CYCLE value = 2.354e-05 2.354e-06 ratio = 10.000
.218800000000000  LIMIT-CYCLE value = 1.548e-05 1.548e-06 ratio = 10.000
.221200000000000  LIMIT-CYCLE value = 1.085e-05 1.085e-06 ratio = 10.000
.223600000000000  LIMIT-CYCLE value = 6.737e-06 6.737e-07 ratio = 10.000
.226000000000000  LIMIT-CYCLE value = 2.413e-05 2.412e-06 ratio = 10.000
.228400000000000  LIMIT-CYCLE value = 1.192e-04 1.192e-05 ratio = 10.000
.230800000000000  LIMIT-CYCLE value = 2.505e-05 2.505e-06 ratio = 10.000
.233200000000000  LIMIT-CYCLE value = 1.082e-05 1.082e-06 ratio = 10.000
.235600000000000  LIMIT-CYCLE value = 7.616e-06 7.616e-07 ratio = 10.000
.238000000000000  LIMIT-CYCLE value = 6.642e-06 6.642e-07 ratio = 10.000
.240400000000000  LIMIT-CYCLE value = 6.029e-06 6.029e-07 ratio = 10.000
.242800000000000  LIMIT-CYCLE value = 5.810e-06 5.810e-07 ratio = 10.000
.245200000000000  LIMIT-CYCLE value = 5.887e-06 5.887e-07 ratio = 10.000
.247600000000000  LIMIT-CYCLE value = 5.689e-06 5.689e-07 ratio = 10.000
real 71.83
user 128.30
sys 4.48
```
The data file can be plotted using plotChaos.py.
```
./plotChaos.py </tmp/results
```

## cns script - Clean Numerical Simulation

This is a relatively new approach to dealing with the global error of ODE simulations, described in detail here: https://arxiv.org/abs/1109.0130.
As a rough guide to the accuracy of a solution, it can be compared with a "better" solution (one made with "better" solver parameters), and discarded at the point where the solutions diverge.
The two simulations are run in parallel processes, but obviously the "better" solution takes longer.
There are three alternative strategies for creating the "better" integrator:

cns shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
step | The step size is reduced by a quarter (suitable for RK4) comparisons
order | The Taylor Series order is increased by two
both | The order is increased by one, and the step size by one half

#### Example output - 30 time units
```
$ ./cns both ./tsm-lorenz-dbg NA NA 10 .01 35000 -15.8 -17.48 35.64 10 28 8 3
Better: ./tsm-lorenz-dbg NA NA 11 .005000 70000 -15.8 -17.48 35.64 10 28 8 3
Divergence: ['./divergence.py', '/tmp/dataA', '/tmp/dataB', '3', '0.000000000001', '0.000000001', '0.000001', '0.001', '1.0']
Threshold: 1.0e-12, t: 0.020
Threshold: 1.0e-09, t: 3.700
Threshold: 1.0e-06, t: 9.160
Threshold: 1.0e-03, t: 20.040
Threshold: 1.0e+00, t: 24.940
```
(matplotlib plot not shown!)

60 time units
```
$ ./cns both ./tsm-lorenz-dbg NA NA 20 .01 65000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

150 time units
```
$ ./cns both ./tsm-lorenz-dbg NA NA 50 .005 150000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

## ic script - Sensitivity to variation in initial conditions
This script is used to generate deviation data for chaos scanning, but the data can also be plotted in real time using matplotlib.
As well as the trajectory specified in the command arguments, six others are created and evolved; each one is the centre of the face of a cube around the original value

ic shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
separation | Initial separation between "original" trajectory and the additional ones
plot / noplot | Allows non-interactive use without plotting the data (e.g. by the chaos scanner).

The simulation is run seven times in parallel processes, the original along with each perturbed x, y, z.
```
$ ./ic .001 plot ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
oo ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
x+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.799 -17.48 35.64 10 28 8 3
x- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.801 -17.48 35.64 10 28 8 3
y+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.479 35.64 10 28 8 3
y- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.481 35.64 10 28 8 3
z+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.641 10 28 8 3
z- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.639 10 28 8 3
```
(3D plot not shown)
