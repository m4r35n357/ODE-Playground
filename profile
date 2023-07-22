#!/bin/sh

args="$0 $*"
echo "args: \033[1;37m$(($# + 1))\033[0;37m, [ \033[0;35m$args\033[0;37m ]" >&2

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
