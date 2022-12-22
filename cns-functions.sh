#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

fileA="/tmp/$USER/dataA"  # results of the "better" simulation
fileB="/tmp/$USER/dataB"  # results of the requested simulation

halfstep () {  # step / 2
    case $1 in
        *tsm-*-std) start="$1 $2 $3";;
                 *) echo "CNS: Bad Executable Name";
                    exit;;
    esac
    step=$(echo "scale=9; $4 / 2;" | /usr/bin/bc)
    steps=$(($5 * 2))
    shift 5
    end="$*"
    echo 'Better:' $start $step $steps $end >&2
    $start $step $steps $end | /usr/bin/sed -n '1~2p' >$fileA &
}
