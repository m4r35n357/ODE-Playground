#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

fileA="/tmp/$USER/dataA"  # results of the "better" simulation
fileB="/tmp/$USER/dataB"  # results of the requested simulation

get_precision () {
    begin="$1 $2"
    case $3 in
        single)    p='23';;
        double)    p='53';;
        extended)  p='63';;
        quadruple) p='113';;
        octuple)   p='237';;
        *)         p="$3";;
    esac
    shift 3
    echo "$begin $p $*"
}

halfstep () {  # step / 2
    start="$1 $2 $3 $4"
    step=$(echo "scale=6; $5 / 2;" | /usr/bin/bc)
    steps=$(($6 * 2))
    shift 6
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~2p' >$fileA &
}
