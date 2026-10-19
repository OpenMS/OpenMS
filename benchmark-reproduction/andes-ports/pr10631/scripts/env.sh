# usage: source env.sh <arm>   (arm = base | pr | fix)
export OPENMS_ARM=$1
export PATH=${OPENMS_PREFIX:-/opt}/openms-$1/bin:$PATH
export LD_LIBRARY_PATH=${OPENMS_PREFIX:-/opt}/openms-$1/lib
export OPENMS_DATA_PATH=${OPENMS_SRC:-/home/user/OpenMS}/share/OpenMS
