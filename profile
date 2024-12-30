#!/bin/sh
#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

. ./base.sh

method=$1
prog=$2
shift 2

echo "\033[33;1mBuilding, please wait . . .\033[0m"
case $method in
     'gcc') make clean && make CCC=prof >/dev/null
           $prog $* >/dev/null
           gprof -p -b $prog gmon.out
           ;;
     'gpt') make clean && make CCC=gpt >/dev/null
           export CPUPROFILE=/tmp/$USER/prof.out
           $prog $* >/dev/null
           google-pprof --text $prog $CPUPROFILE
           ;;
        *) echo "invalid method";;
esac
