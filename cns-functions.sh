#
#  (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

fileA="/tmp/$USER/dataA"  # results of the "better" simulation
fileB="/tmp/$USER/dataB"  # results of the requested simulation

halfstep () {  # step / 2
    start="$1 $2 $3 $4"
    step=$(echo "scale=6; $5 / 2;" | /usr/bin/bc)
    steps=$(($6 * 2))
    shift 6
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~2p' >$fileA &
}

quarterstep () {  # step / 4
    start="$1 $2 $3 $4"
    step=$(echo "scale=6; $5 / 4;" | /usr/bin/bc)
    steps=$(($6 * 4))
    shift 6
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~4p' >$fileA &
}

orderplus1 () {  # order + 1
    start="$1 $2 $3"
    order=$(($4 + 1))
    shift 4
    end="$*"
    echo 'Better:' $start $order $end >&2
    $start $order $end >$fileA &
}

orderstep () {  # order + 1, step / 2
    start="$1 $2 $3"
    order=$(($4 + 1))
    step=$(echo "scale=6; $5 / 2;" | /usr/bin/bc)
    steps=$(($6 * 2))
    shift 6
    end="$*"
    echo 'Better:' $start $order $step $steps $end >&2
    $start $order $step $steps $end | sed -n '1~2p' >$fileA &
}
