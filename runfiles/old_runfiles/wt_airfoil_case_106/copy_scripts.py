copy_names = [
'c105_t21_l15_k16_g2000_n400',
]

import shutil
for name in copy_names: 
    shutil.copyfile('common_runner.py',name + '.py')

# mpirun -n 188 python -m mpi4py c105_t21_l15_k16_g2000_n400.py