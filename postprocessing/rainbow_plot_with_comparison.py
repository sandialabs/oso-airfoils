# This script can be used to add comparison airfoil data to the rainbow plots in the run directories
# Neuralfoil is used presently, but xfoil can be substituted in.  
#
# a 'figurePath' input can be provided to change where the figure is written.

import pathlib
from compare_airfoils import rainbow_plot
path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso = path_to_here.parent
path_to_datfiles = path_to_oso / 'historical_airfoils/mhkf1/'

path_to_data = path_to_here / "cases/cases_111_to_120/case_114/c114_t18_k16_n752_l13_e15__2026_05_13_18-25-5282/population_c114_t18_k16_n752_l13_e15_g500.json"
comp = {"mhkf1-180":str(path_to_datfiles / 'mhkf1-180.dat')}
cor = {"mhkf1-180":'k'}
rainbow_plot(path_to_data, comparison_airfoil=comp, color_override=cor)#, tools = ['xfoil'])

path_to_data = path_to_here / "cases/cases_111_to_120/case_114/c114_t24_k16_n752_l13_e15__2026_05_14_09-51-5547/population_c114_t24_k16_n752_l13_e15_g500.json"
comp = {"mhkf1-240":str(path_to_datfiles / 'mhkf1-240.dat')}
cor = {"mhkf1-240":'k'}
rainbow_plot(path_to_data, comparison_airfoil=comp, color_override=cor)#, tools = ['xfoil'])