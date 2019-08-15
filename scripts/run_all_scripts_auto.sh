#!/bin/bash
module load python/3.6-anaconda-4.4

paramfile=../settings/params_auto_lsst_all_allsimple_biaszmax50.pkl
#srun -n 100 python -u clphiphi_VB_parallel_split.py $paramfile
#srun -n 100 python -u clphiphi_VB_parallel_split2.py $paramfile
srun -n 100 python -u clphiphi_VB_parallel.py $paramfile
#srun -n 25 python -u clpsigal_VB_parallel1.py $paramfile
#srun -n 25 python -u clpsigal_VB_parallel2.py $paramfile
#srun -n 50 python -u clphidelta_VB_parallel_M31b.py $paramfile



