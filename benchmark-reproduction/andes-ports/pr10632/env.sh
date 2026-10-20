# paths; BIN selects the OpenMS build (baseline or PR)
export BIN=${BIN:-/opt/openms-pr/bin}
export LD_LIBRARY_PATH=$(dirname $BIN)/lib
export OPENMS_DATA_PATH=${OPENMS_DATA_PATH:-/home/user/OpenMS/share/OpenMS}
export DATA=${DATA:-/home/user/data}
export WORK=${WORK:-/home/user/work}
export TP=${TP:-/home/user/thirdparty/Linux/x86_64}
export RUNS="LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01 LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02 LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02 LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03"
