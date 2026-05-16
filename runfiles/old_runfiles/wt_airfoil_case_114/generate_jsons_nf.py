import copy
import json

# Resulting call is:
# mpirun -n 8 python -m mpi4py common_runner.py run_001.json
# and similar

default_dict = {
    "case_number"                                :  112,
    "tau"                                        :  0.24,
    "N_k"                                        :  16,
    "N_pop"                                      :  752,
    "N_generations"                              :  500,
    # "CL"                                         :  1.4,
    # "Re"                                         :  13.0e+6,
    "tool"                                       :  "neuralfoil",
    "continuation_file"                          :  None,
    "continuation_file_overwrite"                :  False,
    "outfile_leader"                             :  "./",
    # "xfoil_path"                                 :  None,
    # "xfoil_tempfile_path_leader"                 :  "t_",
    "N_tries"                                    :  1,
    "CMc_min"                                    :  None,
    "CMr_min"                                    :  None,
    "cm_alpha_band"                              :  5.0,
    # "TE_gap"                                     :  0.00751,
    # "Ixx_con"                                    :  0.00041096,
    # "Iyy_con"                                    :  0.00561409,
    # "Izz_con"                                    :  0.00602287,
    # "A_con"                                      :  0.13051205,
    # "ler_con"                                    :  0.025,
    "min_radius_location_upper"                  :  None,
    "min_radius_location_lower"                  :  None,
    "min_radius_location_cutoff"                 :  0.08,
    "cl_max_limit_clean"                         :  None,
    "cl_max_limit_rough"                         :  None,
    "target_stall_margin"                        :  4.0,
    "percent_delta_cl_from_roughness_threshold"  :  0.10,
    "percent_LoD_falloff_threshold"              :  0.15,
    "alpha_falloff_offset"                       :  1.0,
    "max_thickness_loc"                          :  0.275,
    "max_thickness_loc_upper"                    :  0.275,
    "max_thickness_loc_lower"                    :  0.275,
    "ler_skew_factor"                            :  1.9,
    # "cone_angle"                                 :  5.0,
    "te_frac"                                    :  0.95,
    "target_cl"                                  :  None,
    "target_alpha"                               :  None,
    "curvature_bound"                            :  -750,
    "ec_cutoff"                                  :  0.9,
    "cp_min_design"                              :  None,
    "cp_min_at_alpha_offset"                     :  None,
    "cp_min_alpha_offset"                        :  None,
    "cp_min_prestall"                            :  None,
    "stall_margin_clean_weighting"               :  1.0e+2,
    "stall_margin_rough_weighting"               :  1.0e+2,
    "lift_margin_clean_weighting"                :  0.5,
    "cl_max_limit_clean_weighting"               :  1.0e+4,
    "cl_max_limit_rough_weighting"               :  1.0e+4,
    "delta_cl_from_roughness_weighting"          :  1.0e+4,
    "LoD_falloff_weighting"                      :  50.0,
    "ixx_weighting"                              :  1.0e+6,
    "iyy_weighting"                              :  1.0e+4,
    "izz_weighting"                              :  1.0e+4,
    "a_weighting"                                :  1.0e+4,
    "leading_edge_radius_upper_weighting"        :  1.0e+3,
    "leading_edge_radius_lower_weighting"        :  1.0e+3,
    "min_radius_location_upper_weighting"        :  1.0e+4,
    "min_radius_location_lower_weighting"        :  1.0e+4,
    "max_thickness_weighting"                    :  1.0e+4,
    "max_thickness_upper_weighting"              :  1.0e+4,
    "max_thickness_lower_weighting"              :  5.0e+4,
    "radii_skew_weighting"                       :  1.0e+3,
    "curvature_weighting"                        :  100,
    "lower_surface_curvature_weighting"          :  1.0e+2,
    "te_cone_violation_weighting"                :  1.0e+5,
    "CL_target_weighting"                        :  1.0e+4,
    "clean_moment_weighting"                     :  1.0e+4,
    "rough_moment_weighting"                     :  1.0e+4,
    "cp_min_design_weighting"                    :  1.0e+4,
    "cp_min_at_alpha_offset_weighting"           :  1.0e+4,
    "cp_min_prestall_weighting"                  :  1.0e+4,
    "infeasibility_penalty"                      :  1.0e+4,
    "N_crit_clean"                               :  9.0,
    "xtp_u_clean"                                :  1.0,
    "xtp_l_clean"                                :  1.0,
    "alpha_min_clean"                            :  0,
    "alpha_max_clean"                            :  30,
    "alpha_step_clean"                           :  1,
    "N_crit_rough"                               :  3.0,
    "xtp_u_rough"                                :  0.05,
    "xtp_l_rough"                                :  0.05,
    "alpha_min_rough"                            :  0,
    "alpha_max_rough"                            :  20,
    "alpha_step_rough"                           :  1,
    "xfoil_timelimit"                            :  15,
    "neuralfoil_model"                           :  "xxlarge",
    "N_points_moi"                               :  20,
    "N_crossovers"                               :  3,
    "probability_of_mutation"                    :  0.3,
    "maximum_parent_fraction"                    :  0.7,
    "front1_cap_fraction"                        :  0.5,
    "N_mutations"                                :  4,

}


grid_paths = [["/gpfs/ahsieh/tempfiles/xfoil",     "/gpfs/ahsieh/tempfiles/tmp_"    ],
              ["/pscratch/ahsieh/tempfiles/xfoil", "/pscratch/ahsieh/tempfiles/tmp_"],
              ["/tscratch/ahsieh/tempfiles/xfoil", "/tscratch/ahsieh/tempfiles/tmp_"]]

tau_data_dict = {
    0.15: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00196, "Ixx_con": 0.00011000, "Iyy_con": 0.00397999, "Izz_con": 0.00408809, "A_con": 0.08700496, "ler_con": 0.007, "cone_angle": 10.0},
    0.18: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00230, "Ixx_con": 0.00017438, "Iyy_con": 0.00436351, "Izz_con": 0.00454606, "A_con": 0.09995900, "ler_con": 0.008, "cone_angle": 10.0},
    0.21: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00262, "Ixx_con": 0.00027518, "Iyy_con": 0.00493714, "Izz_con": 0.00521632, "A_con": 0.11477620, "ler_con": 0.010, "cone_angle":  5.0},
    0.24: {"CL":1.4, "Re": 13e6, "TE_gap": 0.00751, "Ixx_con": 0.00041096, "Iyy_con": 0.00561409, "Izz_con": 0.00602287, "A_con": 0.13051205, "ler_con": 0.025, "cone_angle":  5.0},
    0.27: {"CL":1.3, "Re": 16e6, "TE_gap": 0.01012, "Ixx_con": 0.00058321, "Iyy_con": 0.00633417, "Izz_con": 0.00691323, "A_con": 0.14660942, "ler_con": 0.030, "cone_angle":  5.0},
    0.30: {"CL":1.2, "Re": 18e6, "TE_gap": 0.01140, "Ixx_con": 0.00079640, "Iyy_con": 0.00706380, "Izz_con": 0.00785849, "A_con": 0.16289864, "ler_con": 0.040, "cone_angle":  0.0},
    0.33: {"CL":1.2, "Re": 16e6, "TE_gap": 0.01140, "Ixx_con": 0.00105795, "Iyy_con": 0.00779600, "Izz_con": 0.00885328, "A_con": 0.17959744, "ler_con": 0.060, "cone_angle":  0.0},
    0.36: {"CL":1.2, "Re": 13e6, "TE_gap": 0.01140, "Ixx_con": 0.00137822, "Iyy_con": 0.00855043, "Izz_con": 0.00991577, "A_con": 0.19731100, "ler_con": 0.080, "cone_angle":  0.0},
}

# CL_vals = [0.8, 1.0, 1.2, 1.4]
CL_vals = [1.0, 1.5]
Re_vals = [5.0, 15.0]
# CL_vals = [1.4]
fnames = []


runstring = ''

counter = 1
for tau in [0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
# for tau in [0.24]:
    write_dict = copy.deepcopy(default_dict)
    write_dict['tau'] = tau

    tau_data = tau_data_dict[tau]
    for ky, vl in tau_data.items():
        write_dict[ky] = vl

    # Need to give the thicker airfoils a little boost in the rough condition
    if tau == 0.36:
        write_dict['target_cl'] = tau_data_dict[tau]['CL'] + 0.3
        write_dict['target_alpha'] = 10.0

    fname = "run_nf_%s.json"%(str(counter).zfill(3))
    fnames.append(fname)
    runstring += fname + ', tau: %.2f\n' % (tau)
    with open(fname, "w") as f:
        rem = counter%3
        write_dict['xfoil_path'] = None #grid_paths[rem][0]
        write_dict['xfoil_tempfile_path_leader'] = "t_" #grid_paths[rem][1]
        json.dump(write_dict, f, indent=4)
    counter += 1

    # # dump to json
    # with open(f"run_%s.json"%(str(counter).zfill(3)), "w") as f:
    #     json.dump(write_dict, f, indent=4)
    # counter += 1

proc_count = None
if write_dict['tool'] == 'xfoil':
    proc_count = 188
else:
    proc_count = 8 #96

shell_text = ''
for fname in fnames:
    shell_text += 'mpirun -n %d python -m mpi4py common_runner.py %s\n' % (proc_count, fname)

f = open("all_run_nf.sh", "w")
f.write(shell_text)
f.close()



f = open("runcases.txt", "w")
f.write(runstring)
f.close()

print(runstring)