copy_names = [
'c101_t21_l15_r110_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r111_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r112_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r113_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r114_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r115_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r116_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r117_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r118_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r119_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r120_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r121_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r122_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r123_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r124_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r125_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r126_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r127_k16_g2000_n400_x5_s3_m_p',
'c101_t21_l15_r128_k16_g2000_n400_x5_s1_m_p',
'c101_t21_l15_r129_k16_g2000_n400_x5_s2_m_p',
'c101_t21_l15_r130_k16_g2000_n400_x5_s3_m_p',
]

import shutil
for name in copy_names: 
    shutil.copyfile('common_runner.py',name + '.py')

# mpirun -n 188 python -m mpi4py c101_t21_l15_r122_k16_g2000_n400_x5_s1_m_p.py