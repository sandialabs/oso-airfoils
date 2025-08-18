copy_names = [
'c67_t15_l15_r124_k16_n400',
'c67_t18_l15_r123_k16_n400',
'c67_t21_l10_r122_k16_n400',
'c67_t21_l11_r122_k16_n400',
'c67_t21_l12_r122_k16_n400',
'c67_t21_l13_r122_k16_n400',
'c67_t21_l14_r122_k16_n400',
'c67_t21_l15_r115_k16_n400',
'c67_t21_l15_r116_k16_n400',
'c67_t21_l15_r117_k16_n400',
'c67_t21_l15_r118_k16_n400',
'c67_t21_l15_r119_k16_n400',
'c67_t21_l15_r120_k16_n400',
'c67_t21_l15_r121_k16_n400',
'c67_t21_l15_r122_k16_n400',
'c67_t24_l14_r118_k16_n400',
'c67_t27_l13_r111_k16_n400',
'c67_t30_l12_r98_k16_n400',
'c67_t33_l12_r90_k16_n400',
'c67_t36_l12_r78_k16_n400',
'c67_t24_l14_r118_k16_n376',
]

import shutil
for name in copy_names: 
    shutil.copyfile('common_runner.py',name + '.py')

# mpirun -n 188 python -m mpi4py c64_t24_l14_r118_k16_n388.py