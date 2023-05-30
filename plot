#!/bin/sh
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

args="$0 $*"
echo "args: \033[1;37m$(($# + 1))\033[0;37m, [ \033[0;35m$args\033[0;37m ]" >&2

user_data="/tmp/$USER/data"

plot () {
 gnuplot -p << EOF
set terminal wxt background rgb "grey85"
set xyplane 0
set view 54.73561,135
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'
splot '$user_data' with lines
EOF
}

$* >$user_data
plot
