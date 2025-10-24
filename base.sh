#
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

GRY='\033[1;30m'
RED='\033[1;31m'
GRN='\033[1;32m'
YLW='\033[1;33m'
BLU='\033[1;34m'
MGT='\033[0;35m'
CYN='\033[0;36m'
WHT='\033[1;37m'
NRM='\033[0m'

args="$0 $*"
echo "${GRY}args ${NRM}$(($# + 1))${GRY}, argv [ ${MGT}${args}${GRY} ]${NRM}" >&2

user_dir="/tmp/$USER"
[ -d $user_dir ] || mkdir $user_dir
user_data="$user_dir/data"
