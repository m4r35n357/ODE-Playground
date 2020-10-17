#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

quarterstep () {  # step / 4
    start="$1 $2 $3 $4"
    step=$(echo "scale=6; $5 / 4;" | bc)
    steps=$(($6 * 4))
    shift 6
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~4p' >$fileA &
}

eightthstep () {  # step / 8
    start="$1 $2 $3 $4"
    step=$(echo "scale=6; $5 / 8;" | bc)
    steps=$(($6 * 8))
    shift 6
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | sed -n '1~8p' >$fileA &
}

orderplus2 () {  # order + 2
    start="$1 $2 $3"
    order=$(($4 + 2))
    shift 4
    end="$*"
    echo 'Better:' $start $order $end >&2
    $start $order $end >$fileA &
}

orderstep () {  # order + 1, step / 2
    start="$1 $2 $3"
    order=$(($4 + 1))
    step=$(echo "scale=6; $5 / 2;" | bc)
    steps=$(($6 * 2))
    shift 6
    end="$*"
    echo 'Better:' $start $order $step $steps $end >&2
    $start $order $step $steps $end | sed -n '1~2p' >$fileA &
}

orderstep2 () {  # order + 2, step / 4
    start="$1 $2 $3"
    order=$(($4 + 2))
    step=$(echo "scale=6; $5 / 4;" | bc)
    steps=$(($6 * 4))
    shift 6
    end="$*"
    echo 'Better:' $start $order $step $steps $end >&2
    $start $order $step $steps $end | sed -n '1~4p' >$fileA &
}
