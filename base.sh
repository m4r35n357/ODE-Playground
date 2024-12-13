#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

GRY='\e[1;30m'
RED='\e[1;31m'
GRN='\e[1;32m'
YLW='\e[1;33m'
BLU='\e[1;34m'
MGT='\e[0;35m'
CYN='\e[0;36m'
WHT='\e[1;37m'
NRM='\e[0m'

args="$0 $*"
echo "${GRY}args ${NRM}$(($# + 1))${GRY}, argv [ ${MGT}${args}${GRY} ]${NRM}" >&2

user_dir="/tmp/$USER"
[ -d $user_dir ] || mkdir $user_dir
user_data="$user_dir/data"
