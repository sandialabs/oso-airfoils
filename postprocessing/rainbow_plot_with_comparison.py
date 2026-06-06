# This script can be used to add comparison airfoil data to the rainbow plots in the run directories
# Neuralfoil is used presently, but xfoil can be substituted in.  
#
# a 'figurePath' input can be provided to change where the figure is written.

import pathlib
from compare_airfoils import rainbow_plot
path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso = path_to_here.parent
path_to_datfiles = path_to_oso / 'historical_airfoils/mhkf1/'

# # path_to_data = path_to_here / "cases/cases_111_to_120/case_115/c115_t18_k16_n752_l13_e15__2026_05_19_23-59-1193/population_c115_t18_k16_n752_l13_e15_g065.json"
path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t18_k16_n752_l13_e15__2026_06_04_04-56-3876/population_c115_t18_k16_n752_l13_e15_g500.json"
comp = {"mhkf1-180":str(path_to_datfiles / 'mhkf1-180.dat')}
cor = {"mhkf1-180":'k'}
rainbow_plot(path_to_data, comparison_airfoil=comp, color_override=cor)#, tools = ['xfoil'])

# path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t18_k16_n752_l13_e15__2026_05_21_00-51-0784/population_c115_t18_k16_n752_l13_e15_g500.json"
# comp = {"mhkf1-180":str(path_to_datfiles / 'mhkf1-180.dat')}
# cor = {"mhkf1-180":'k'}
# rainbow_plot(path_to_data, comparison_airfoil=comp, color_override=cor)#, tools = ['xfoil'])

# path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t21_k16_n752_l13_e15__2026_05_20_08-06-4000/population_c115_t21_k16_n752_l13_e15_g500.json"
# rainbow_plot(path_to_data, comparison_airfoil={}, color_override={})#, tools = ['xfoil'])

# path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t21_k16_n752_l13_e15__2026_05_21_08-38-0187/population_c115_t21_k16_n752_l13_e15_g500.json"
# rainbow_plot(path_to_data, comparison_airfoil={}, color_override={})#, tools = ['xfoil'])

# # path_to_data = path_to_here / "cases/cases_111_to_120/case_114/c114_t24_k16_n752_l13_e15__2026_05_14_09-51-5547/population_c114_t24_k16_n752_l13_e15_g500.json"
# path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t24_k16_n752_l13_e15__2026_05_20_17-04-3892/population_c115_t24_k16_n752_l13_e15_g500.json"
# comp = {"mhkf1-240":str(path_to_datfiles / 'mhkf1-240.dat')}
# cor = {"mhkf1-240":'k'}
# rainbow_plot(path_to_data, comparison_airfoil=comp, color_override=cor)#, tools = ['xfoil'])


# path_to_data = "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_115/c115_t33_k16_n752_l12_e15__2026_05_22_15-47-1615/population_c115_t33_k16_n752_l12_e15_g500.json"
# rainbow_plot(path_to_data, comparison_airfoil={}, color_override={})#, tools = ['xfoil'])