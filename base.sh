
args="$0 $*"
echo "args \033[1;37m$(($# + 1))\033[0;37m [ \033[0;35m$args\033[0;37m ]" >&2

user_dir="/tmp/$USER"
[ -d $user_dir ] || mkdir $user_dir
user_data="$user_dir/data"

