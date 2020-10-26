#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

fileA='/tmp/dataA'  # results of the "better" simulation
fileB='/tmp/dataB'  # results of the requested simulation

quarterstep () {  # step / 4
    start="$1 $2 $3"
    step=$(echo "scale=6; $4 / 4;" | /usr/bin/bc)
    steps=$(($5 * 4))
    shift 5
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~4p' >$fileA &
}

eightthstep () {  # step / 8
    start="$1 $2 $3"
    step=$(echo "scale=6; $4 / 8;" | /usr/bin/bc)
    steps=$(($5 * 8))
    shift 5
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~8p' >$fileA &
}

orderplus2 () {  # order + 2
    start="$1 $2"
    order=$(($3 + 2))
    shift 3
    end="$*"
    echo 'Better:' $start $order $end >&2
    $start $order $end >$fileA &
}

orderstep () {  # order + 1, step / 2
    start="$1 $2"
    order=$(($3 + 1))
    step=$(echo "scale=6; $4 / 2;" | /usr/bin/bc)
    steps=$(($5 * 2))
    shift 5
    end="$*"
    echo 'Better:' $start $order $step $steps $end >&2
    $start $order $step $steps $end | sed -n '1~2p' >$fileA &
}

orderstep2 () {  # order + 2, step / 4
    start="$1 $2"
    order=$(($3 + 2))
    step=$(echo "scale=6; $4 / 4;" | /usr/bin/bc)
    steps=$(($5 * 4))
    shift 5
    end="$*"
    echo 'Better:' $start $order $step $steps $end >&2
    $start $order $step $steps $end | sed -n '1~4p' >$fileA &
}
