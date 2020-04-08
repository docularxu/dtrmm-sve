#!/bin/bash
echo `pwd`
#For Mac OS only
export DYLD_LIBRARY_PATH=/opt/intel/lib:/opt/intel/mkl/lib

#Single Thread
export KMP_AFFINITY=compact  #Rule to bind core to thread for OMP thread with Intel compiler for parallel version
export OMP_NUM_THREADS=1     #Set OMP number of threads for parallel version
export BLISLAB_IC_NT=1       #Set BLISLAB number of threads for parallel version
k_start=16
k_end=10240
k_blocksize=16
echo "result=["
echo -e "%m\t%n\t%k\t%MY_GFLOPS\t%REF_GFLOPS"
for (( k=k_start; k<=k_end; k+=k_blocksize ))
do
    ./test_bl_dtrmm.x     $k $k
done
echo "];"

