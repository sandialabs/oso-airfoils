copy_names = [
'c73_t21_l15_r122_k16_g1000_n400_x5_s1',
'c73_t24_l14_r118_k16_g1000_n400_x5_s2',
'c73_t27_l13_r111_k16_g1000_n400_x5_s3',
]

import shutil
for name in copy_names: 
    shutil.copyfile('common_runner.py',name + '.py')

# mpirun -n 188 python -m mpi4py c73_t21_l15_r122_k16_g1000_n376_x1_m_p.py