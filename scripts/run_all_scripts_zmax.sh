#!/bin/bash
module load python/3.6-anaconda-4.4
#only one run per zmax!
#paramfile=../settings/params_lsstall_cmblens_simple_bias_zmax10886.pkl
#srun -n 100 python -u clphiphi_VB_parallel_split.py $paramfile
#srun -n 100 python -u clphiphi_VB_parallel_split2.py $paramfile
#srun -n 100 python -u clphiphi_VB_parallel.py $paramfile
#srun -n 100 python -u clpsiphi_VB_parallel.py $paramfile
#srun -n 100 python -u clphidelta_VB_parallel_MB2.py $paramfile
paramfile=../settings/params_gaussgal_z10_sigma1_deltalens_z7_simple_bias_zmax15.pkl
srun -n 100 python -u clphiphi_VB_parallel_split.py $paramfile
srun -n 100 python -u clphiphi_VB_parallel_split2.py $paramfile
srun -n 100 python -u clphiphi_VB_parallel.py $paramfile
srun -n 100 python -u clpsiphi_VB_parallel.py $paramfile
srun -n 100 python -u clphidelta_VB_parallel_MB2.py $paramfile

