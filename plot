#!/bin/sh
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

. ./base.sh

keep=$1
shift

plot3d () {
/usr/bin/gnuplot << EOF
set terminal qt
set title '$args'
set xyplane 0
set view 54.73561,135
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'
splot '$user_data' with lines
pause mouse close
EOF
}

$* | tail -n $keep >$user_data

plot3d &
