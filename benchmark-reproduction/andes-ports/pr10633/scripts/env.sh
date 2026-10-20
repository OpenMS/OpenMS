# Common settings. Source this file.
export BENCH=/home/user/bench
export DATA=$BENCH/data
export DATA_PB=$BENCH/data_pb
export RUNS="LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01 LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02 LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03"
export OMP_WAIT_POLICY=passive
# use the OpenMS installation in $1 (bin/ lib/ share/)
use_openms() {
  export OPENMS_ROOT=$1
  export PATH=$1/bin:/usr/bin:/bin
  export LD_LIBRARY_PATH=$1/lib
  export OPENMS_DATA_PATH=$1/share/OpenMS
}
