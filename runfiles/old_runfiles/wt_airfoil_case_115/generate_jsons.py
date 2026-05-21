import copy
import json

# Resulting call is:
# mpirun -n 8 python -m mpi4py common_runner.py run_001.json
# and similar

default_dict = {
    "case_number"                                :  115,
    # "tau"                                        :  0.24,
    "N_k"                                        :  16,
    "N_pop"                                      :  752,
    "N_generations"                              :  500,
    # "CL"                                         :  1.4,
    # "Re"                                         :  13.0e+6,
    "tool"                                       :  "xfoil",
    "continuation_file"                          :  None,
    "continuation_file_overwrite"                :  False,
    "outfile_leader"                             :  "./",
    # "xfoil_path"                                 :  None,
    # "xfoil_tempfile_path_leader"                 :  "t_",
    "N_tries"                                    :  3,
    "CMc_min"                                    :  None,
    "CMr_min"                                    :  None,
    "cm_alpha_band"                              :  5.0,
    # "TE_gap"                                     :  0.00751,
    # "Ixx_con"                                    :  0.00041096,
    # "Iyy_con"                                    :  0.00561409,
    # "Izz_con"                                    :  0.00602287,
    # "A_con"                                      :  0.13051205,
    # "ler_con"                                    :  0.025,
    "min_radius_location_upper"                  :  None, # applies default
    "min_radius_location_lower"                  :  None, # applies default
    "min_radius_location_cutoff"                 :  0.08,
    "cl_max_limit_clean"                         :  None,
    "cl_max_limit_rough"                         :  None,
    "target_stall_margin"                        :  4.0,
    # "percent_delta_cl_from_roughness_threshold"  :  0.10,
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
    # "cp_min_design"                              :  None,
    "cp_min_at_alpha_offset"                     :  None,
    "cp_min_alpha_offset"                        :  None,
    # "cp_min_prestall"                            :  None,
    "toothpick_height"                           :  0.00, 
    "toothpick_location"                         :  0.80,
    "stall_margin_clean_weighting"               :  1.0e+2,
    "stall_margin_rough_weighting"               :  1.0e+2,
    "lift_margin_clean_weighting"                :  0.5,
    "cl_max_limit_clean_weighting"               :  None,
    "cl_max_limit_rough_weighting"               :  None,
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
    "CL_target_weighting"                        :  None,
    "clean_moment_weighting"                     :  None,
    "rough_moment_weighting"                     :  None,
    # "cp_min_design_weighting"                    :  1.0e+4,
    "cp_min_at_alpha_offset_weighting"           :  None,
    # "cp_min_prestall_weighting"                  :  1.0e+4,
    "toothpick_weighting"                        :  None,
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
    "maximum_parent_fraction"                    :  0.3,
    "front1_cap_fraction"                        :  0.5,
    "N_mutations"                                :  4,

}


grid_paths = [["/gpfs/ahsieh/tempfiles/xfoil",     "/gpfs/ahsieh/tempfiles/tmp_"    ],
              ["/pscratch/ahsieh/tempfiles/xfoil", "/pscratch/ahsieh/tempfiles/tmp_"],
              ["/tscratch/ahsieh/tempfiles/xfoil", "/tscratch/ahsieh/tempfiles/tmp_"]]


wt_cl_vals = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
tau_data_dict_wt = {
    0.15: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00196, "Ixx_con": 0.00011000, "Iyy_con": 0.00397999, "Izz_con": 0.00408809, "A_con": 0.08700496, "ler_con_upper": 0.007, "ler_con_lower": 0.007, "cone_angle": 10.0},
    0.18: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00230, "Ixx_con": 0.00017438, "Iyy_con": 0.00436351, "Izz_con": 0.00454606, "A_con": 0.09995900, "ler_con_upper": 0.008, "ler_con_lower": 0.008, "cone_angle": 10.0},
    0.21: {"CL":1.5, "Re": 12e6, "TE_gap": 0.00262, "Ixx_con": 0.00027518, "Iyy_con": 0.00493714, "Izz_con": 0.00521632, "A_con": 0.11477620, "ler_con_upper": 0.010, "ler_con_lower": 0.010, "cone_angle":  5.0},
    0.24: {"CL":1.4, "Re": 13e6, "TE_gap": 0.00751, "Ixx_con": 0.00041096, "Iyy_con": 0.00561409, "Izz_con": 0.00602287, "A_con": 0.13051205, "ler_con_upper": 0.025, "ler_con_lower": 0.025, "cone_angle":  5.0},
    0.27: {"CL":1.3, "Re": 16e6, "TE_gap": 0.01012, "Ixx_con": 0.00058321, "Iyy_con": 0.00633417, "Izz_con": 0.00691323, "A_con": 0.14660942, "ler_con_upper": 0.030, "ler_con_lower": 0.030, "cone_angle":  5.0},
    0.30: {"CL":1.2, "Re": 18e6, "TE_gap": 0.01140, "Ixx_con": 0.00079640, "Iyy_con": 0.00706380, "Izz_con": 0.00785849, "A_con": 0.16289864, "ler_con_upper": 0.040, "ler_con_lower": 0.040, "cone_angle":  0.0},
    0.33: {"CL":1.2, "Re": 16e6, "TE_gap": 0.01140, "Ixx_con": 0.00105795, "Iyy_con": 0.00779600, "Izz_con": 0.00885328, "A_con": 0.17959744, "ler_con_upper": 0.060, "ler_con_lower": 0.060, "cone_angle":  0.0},
    0.36: {"CL":1.2, "Re": 13e6, "TE_gap": 0.01140, "Ixx_con": 0.00137822, "Iyy_con": 0.00855043, "Izz_con": 0.00991577, "A_con": 0.19731100, "ler_con_upper": 0.080, "ler_con_lower": 0.080, "cone_angle":  0.0},
    'all': {
        "cp_min_design"                              :  None,
        "cp_min_prestall"                            :  None,
        "cp_min_design_weighting"                    :  None,
        "cp_min_prestall_weighting"                  :  None,
        "percent_delta_cl_from_roughness_threshold"  :  0.10,
    }
}

# baseline
tau_data_dict_ht_1 = {
    0.18: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00230, "Ixx_con": 0.00017438, "Iyy_con": 0.00436351, "Izz_con": 0.00454606, "A_con": 0.09995900, "ler_con_upper": 0.008, "ler_con_lower": 0.008, "cone_angle": 10.0},
    0.21: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00262, "Ixx_con": 0.00027518, "Iyy_con": 0.00493714, "Izz_con": 0.00521632, "A_con": 0.11477620, "ler_con_upper": 0.010, "ler_con_lower": 0.010, "cone_angle":  5.0},
    0.24: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00751, "Ixx_con": 0.00041096, "Iyy_con": 0.00561409, "Izz_con": 0.00602287, "A_con": 0.13051205, "ler_con_upper": 0.025, "ler_con_lower": 0.025, "cone_angle":  5.0},
    0.27: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.01012, "Ixx_con": 0.00058321, "Iyy_con": 0.00633417, "Izz_con": 0.00691323, "A_con": 0.14660942, "ler_con_upper": 0.030, "ler_con_lower": 0.030, "cone_angle":  5.0},
    0.30: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00079640, "Iyy_con": 0.00706380, "Izz_con": 0.00785849, "A_con": 0.16289864, "ler_con_upper": 0.040, "ler_con_lower": 0.040, "cone_angle":  0.0},
    0.33: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00105795, "Iyy_con": 0.00779600, "Izz_con": 0.00885328, "A_con": 0.17959744, "ler_con_upper": 0.060, "ler_con_lower": 0.060, "cone_angle":  0.0},
    0.36: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00137822, "Iyy_con": 0.00855043, "Izz_con": 0.00991577, "A_con": 0.19731100, "ler_con_upper": 0.080, "ler_con_lower": 0.080, "cone_angle":  0.0},
    'all': {
        "cp_min_design"                              :  -2.0,
        "cp_min_prestall"                            :  -7.0,
        "cp_min_design_weighting"                    :  1.0e+4,
        "cp_min_prestall_weighting"                  :  1.0e+4,
        "percent_delta_cl_from_roughness_threshold"  :  0.10,
    }
}

# radii that match the existing mhkf1 hydrofoils
tau_data_dict_ht_2 = {
    0.18: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00230, "Ixx_con": 0.00017438, "Iyy_con": 0.00436351, "Izz_con": 0.00454606, "A_con": 0.09995900, "ler_con_upper": 0.0263304312898284, "ler_con_lower": 0.0123170677406472, "cone_angle": 10.0},
    0.21: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00262, "Ixx_con": 0.00027518, "Iyy_con": 0.00493714, "Izz_con": 0.00521632, "A_con": 0.11477620, "ler_con_upper": 0.0272062376903069, "ler_con_lower": 0.0208389424521003, "cone_angle":  5.0},
    0.24: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00751, "Ixx_con": 0.00041096, "Iyy_con": 0.00561409, "Izz_con": 0.00602287, "A_con": 0.13051205, "ler_con_upper": 0.0280820440907854, "ler_con_lower": 0.0293608171635534, "cone_angle":  5.0},
    0.27: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.01012, "Ixx_con": 0.00058321, "Iyy_con": 0.00633417, "Izz_con": 0.00691323, "A_con": 0.14660942, "ler_con_upper": 0.0289578504912638, "ler_con_lower": 0.0378826918750065, "cone_angle":  5.0},
    0.30: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00079640, "Iyy_con": 0.00706380, "Izz_con": 0.00785849, "A_con": 0.16289864, "ler_con_upper": 0.0298336568917423, "ler_con_lower": 0.0464045665864596, "cone_angle":  0.0},
    0.33: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00105795, "Iyy_con": 0.00779600, "Izz_con": 0.00885328, "A_con": 0.17959744, "ler_con_upper": 0.0307094632922208, "ler_con_lower": 0.0549264412979128, "cone_angle":  0.0},
    0.36: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00137822, "Iyy_con": 0.00855043, "Izz_con": 0.00991577, "A_con": 0.19731100, "ler_con_upper": 0.0315852696926993, "ler_con_lower": 0.0634483160093659, "cone_angle":  0.0},
    'all': {
        "cp_min_design"                              :  -2.0,
        "cp_min_prestall"                            :  -7.0,
        "cp_min_design_weighting"                    :  1.0e+4,
        "cp_min_prestall_weighting"                  :  1.0e+4,
        "percent_delta_cl_from_roughness_threshold"  :  0.10,
    }
}

# tighten the requirements on cl loss
tau_data_dict_ht_3 = {
    0.18: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00230, "Ixx_con": 0.00017438, "Iyy_con": 0.00436351, "Izz_con": 0.00454606, "A_con": 0.09995900, "ler_con_upper": 0.008, "ler_con_lower": 0.008, "cone_angle": 10.0},
    0.21: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00262, "Ixx_con": 0.00027518, "Iyy_con": 0.00493714, "Izz_con": 0.00521632, "A_con": 0.11477620, "ler_con_upper": 0.010, "ler_con_lower": 0.010, "cone_angle":  5.0},
    0.24: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.00751, "Ixx_con": 0.00041096, "Iyy_con": 0.00561409, "Izz_con": 0.00602287, "A_con": 0.13051205, "ler_con_upper": 0.025, "ler_con_lower": 0.025, "cone_angle":  5.0},
    0.27: {"CL":1.3, "Re": 1.5e6, "TE_gap": 0.01012, "Ixx_con": 0.00058321, "Iyy_con": 0.00633417, "Izz_con": 0.00691323, "A_con": 0.14660942, "ler_con_upper": 0.030, "ler_con_lower": 0.030, "cone_angle":  5.0},
    0.30: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00079640, "Iyy_con": 0.00706380, "Izz_con": 0.00785849, "A_con": 0.16289864, "ler_con_upper": 0.040, "ler_con_lower": 0.040, "cone_angle":  0.0},
    0.33: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00105795, "Iyy_con": 0.00779600, "Izz_con": 0.00885328, "A_con": 0.17959744, "ler_con_upper": 0.060, "ler_con_lower": 0.060, "cone_angle":  0.0},
    0.36: {"CL":1.2, "Re": 1.5e6, "TE_gap": 0.01140, "Ixx_con": 0.00137822, "Iyy_con": 0.00855043, "Izz_con": 0.00991577, "A_con": 0.19731100, "ler_con_upper": 0.080, "ler_con_lower": 0.080, "cone_angle":  0.0},
    'all': {
        "cp_min_design"                              :  -2.0,
        "cp_min_prestall"                            :  -7.0,
        "cp_min_design_weighting"                    :  1.0e+4,
        "cp_min_prestall_weighting"                  :  1.0e+4,
        "percent_delta_cl_from_roughness_threshold"  :  0.05,
    }
}

ctr = 1
fnames = []
runstr = "mpirun -n %d python -m mpi4py common_runner.py %s"
runtext = ""
xfoil_procs = 188
nfoil_procs = 96

for i, dtdct in enumerate([tau_data_dict_ht_1, tau_data_dict_ht_2, tau_data_dict_ht_3], start=1):
    for tau in [0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
        write_dict = copy.deepcopy(default_dict)
        for ky, vl in dtdct['all'].items():
            write_dict[ky] = vl
        spec_dict = copy.deepcopy(dtdct[tau])
        for ky, vl in spec_dict.items():
            write_dict[ky] = vl
        
        write_dict['tau'] = tau
        write_dict['TE_gap'] = 0.0192 * tau/0.18

        if tau == 0.36:
            write_dict['target_cl'] = dtdct[tau]['CL'] + 0.3
            write_dict['target_alpha'] = 10.0
            write_dict['CL_target_weighting'] = 1.0e+5

        nmst = 'base' if i == 1 else ('rad' if i == 2 else 'tight')
        fname = "run_jsons/run_ht_%s_t%s_c%02d_x_%s.json"%(str(ctr).zfill(3), str(int(tau*100)).zfill(2), int(write_dict['CL']*10), nmst)
        fnames.append(fname)
        runtext += runstr%(xfoil_procs, fname) + "\n"
        with open(fname, "w") as f:
            rem = ctr%3
            write_dict['xfoil_path'] = grid_paths[rem][0]
            write_dict['xfoil_tempfile_path_leader'] = grid_paths[rem][1]
            json.dump(write_dict, f, indent=4)

            ctr += 1

for tau in [0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
    for cl in wt_cl_vals:
        write_dict = copy.deepcopy(default_dict)
        for ky, vl in tau_data_dict_wt['all'].items():
            write_dict[ky] = vl
        spec_dict = copy.deepcopy(tau_data_dict_wt[tau])
        for ky, vl in spec_dict.items():
            write_dict[ky] = vl
        write_dict['CL'] = cl

        write_dict['tau'] = tau

        if tau == 0.36:
            write_dict['target_cl'] = tau_data_dict_wt[tau]['CL'] + 0.3
            write_dict['target_alpha'] = 10.0
            write_dict['CL_target_weighting'] = 1.0e+5

        fname = "run_jsons/run_wt_%s_t%s_c%02d_x.json"%(str(ctr).zfill(3), str(int(tau*100)).zfill(2), int(write_dict['CL']*10))
        fnames.append(fname)
        runtext += runstr%(xfoil_procs, fname) + "\n"
        with open(fname, "w") as f:
            rem = ctr%3
            write_dict['xfoil_path'] = grid_paths[rem][0]
            write_dict['xfoil_tempfile_path_leader'] = grid_paths[rem][1]
            json.dump(write_dict, f, indent=4)

        ctr += 1


for i, dtdct in enumerate([tau_data_dict_ht_1, tau_data_dict_ht_2, tau_data_dict_ht_3], start=1):
    for tau in [0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
        write_dict = copy.deepcopy(default_dict)
        for ky, vl in dtdct['all'].items():
            write_dict[ky] = vl
        spec_dict = copy.deepcopy(dtdct[tau])
        for ky, vl in spec_dict.items():
            write_dict[ky] = vl

        write_dict['tau'] = tau
        write_dict['TE_gap'] = 0.0192 * tau/0.18

        if tau == 0.36:
            write_dict['target_cl'] = dtdct[tau]['CL'] + 0.3
            write_dict['target_alpha'] = 10.0
            write_dict['CL_target_weighting'] = 1.0e+5

        write_dict['tool'] = 'neuralfoil'
        write_dict['N_tries'] = 1

        nmst = 'base' if i == 1 else ('rad' if i == 2 else 'tight')
        fname = "run_jsons/run_ht_%s_t%s_c%02d_n_%s.json"%(str(ctr).zfill(3), str(int(tau*100)).zfill(2), int(write_dict['CL']*10), nmst)
        fnames.append(fname)
        runtext += runstr%(nfoil_procs, fname) + "\n"
        with open(fname, "w") as f:
            rem = ctr%3
            write_dict['xfoil_path'] = grid_paths[rem][0]
            write_dict['xfoil_tempfile_path_leader'] = grid_paths[rem][1]
            json.dump(write_dict, f, indent=4)

            ctr += 1


for fname in fnames:
    print(fname)

# runtext = runtext.replace("run_jsons/", "")
f = open("all_run.sh", "w")
f.write(runtext)
f.close()
