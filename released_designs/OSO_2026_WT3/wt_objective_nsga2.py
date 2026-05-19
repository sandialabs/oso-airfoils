import numpy as np
from kulfan import Kulfan, units
import copy
import math

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()
from xfoil_wrapper_noprint import run as run_xfoil
from neuralfoil_wrapper_noprint import run as run_neuralfoil

from geometry_functions import (
    TE_gap_function,
    cone_angle_function,
    Ixx_function,
    Iyy_function,
    Izz_function,
    area_function,
    ler_function,
    min_radius_location_upper_function,
    min_radius_location_lower_function)


def airfoil_fitness(x):
    N_tries = x['params']['N_tries']
    if N_tries is None:
        N_rtr = 1
    else:
        N_rtr = N_tries

    for i in range(0,N_rtr):
        res = core_fitness_function(x)
        if res[3] > -10:
            return res
    return res

def core_fitness_function(x):
    # ----------------------
    # unpack
    # ----------------------
    pid = x['pid']
    K_upper = x['individual'][0:int(len(x['individual'])/2)]
    K_lower = x['individual'][int(len(x['individual'])/2):]
    # te_gap = x['individual'][-1]
    tau         = x['params']['tau']
    CL_in       = x['params']['CL']
    CMc_in      = x['params']['CMc_min']
    CMr_in      = x['params']['CMr_min']
    Re_in       = x['params']['Re']

    N_reported    = 16
    N_constraints = 31

    # for iiiiii in range (0,1):
    try:

        cl_design = CL_in
        Re = Re_in

        if 'TE_gap' not in x['params'] or x['params']['TE_gap'] is None:
            te_gap = TE_gap_function(tau)
        else:
            te_gap = x['params']['TE_gap']
        
        if 'cone_angle' not in x['params'] or x['params']['cone_angle'] is None:
            cone_angle = cone_angle_function(tau)
            te_frac = 0.95
        else:
            cone_angle = x['params']['cone_angle']
            te_frac    = x['params']['te_frac']

        if 'Ixx_con' not in x['params'] or x['params']['Ixx_con'] is None:
            Ixx_con = Ixx_function(tau)
        else:
            Ixx_con = x['params']['Ixx_con']

        if 'Iyy_con' not in x['params'] or x['params']['Iyy_con'] is None:
            Iyy_con = Iyy_function(tau)
        else:
            Iyy_con = x['params']['Iyy_con']

        if 'Izz_con' not in x['params'] or x['params']['Izz_con'] is None:
            Izz_con = Izz_function(tau)
        else:
            Izz_con = x['params']['Izz_con']

        if 'A_con' not in x['params'] or x['params']['A_con'] is None:
            A_con = area_function(tau)
        else:
            A_con = x['params']['A_con']

        if 'ler_con' not in x['params'] or x['params']['ler_con'] is None:
            ler_con = ler_function(tau)
        else:
            ler_con = x['params']['ler_con']
        
        if 'min_radius_location_upper' not in x['params'] or x['params']['min_radius_location_upper'] is None:
            min_radius_location_upper = min_radius_location_upper_function(tau)
        else:
            min_radius_location_upper = x['params']['min_radius_location_upper']

        if 'min_radius_location_lower' not in x['params'] or x['params']['min_radius_location_lower'] is None:
            min_radius_location_lower = min_radius_location_lower_function(tau)
        else:
            min_radius_location_lower = x['params']['min_radius_location_lower']

        if 'min_radius_location_cutoff' not in x['params'] or x['params']['min_radius_location_cutoff'] is None:
            min_radius_location_cutoff = 0.08
        else:
            min_radius_location_cutoff = x['params']['min_radius_location_cutoff']

        target_stall_margin                       = x['params']['target_stall_margin']
        percent_delta_cl_from_roughness_threshold = x['params']['percent_delta_cl_from_roughness_threshold']
        alpha_falloff_offset                      = x['params']['alpha_falloff_offset']
        percent_LoD_falloff_threshold             = x['params']['percent_LoD_falloff_threshold']
        max_thickness_loc                         = x['params']['max_thickness_loc']
        max_thickness_loc_upper                   = x['params']['max_thickness_loc_upper']
        max_thickness_loc_lower                   = x['params']['max_thickness_loc_lower']
        ler_skew_factor                           = x['params']['ler_skew_factor']
        cl_max_limit_clean                        = x['params']['cl_max_limit_clean']
        cl_max_limit_rough                        = x['params']['cl_max_limit_rough']
        if 'target_cl' not in x['params'] or x['params']['target_cl'] is None:
            target_cl = None
        else:
            target_cl = x['params']['target_cl']
        if 'target_alpha' not in x['params'] or x['params']['target_alpha'] is None:
            target_alpha = None
        else:
            target_alpha = x['params']['target_alpha']
        curvature_bound                           = x['params']['curvature_bound']
        ec_cutoff                                 = x['params']['ec_cutoff']
        cm_alpha_band                             = x['params']['cm_alpha_band']
        cp_min_design                             = x['params']['cp_min_design']
        cp_min_at_alpha_offset                    = x['params']['cp_min_at_alpha_offset']
        cp_min_alpha_offset                       = x['params']['cp_min_alpha_offset']
        cp_min_prestall                           = x['params']['cp_min_prestall']
        toothpick_height                          = x['params']['toothpick_height']
        toothpick_location                        = x['params']['toothpick_location']

        stall_margin_clean_weighting        = x['params']['stall_margin_clean_weighting']
        stall_margin_rough_weighting        = x['params']['stall_margin_rough_weighting']
        lift_margin_clean_weighting         = x['params']['lift_margin_clean_weighting']
        cl_max_limit_clean_weighting        = x['params']['cl_max_limit_clean_weighting']
        cl_max_limit_rough_weighting        = x['params']['cl_max_limit_rough_weighting']
        delta_cl_from_roughness_weighting   = x['params']['delta_cl_from_roughness_weighting']
        LoD_falloff_weighting               = x['params']['LoD_falloff_weighting']
        ixx_weighting                       = x['params']['ixx_weighting']
        iyy_weighting                       = x['params']['iyy_weighting']
        izz_weighting                       = x['params']['izz_weighting']
        a_weighting                         = x['params']['a_weighting']
        leading_edge_radius_upper_weighting = x['params']['leading_edge_radius_upper_weighting']
        leading_edge_radius_lower_weighting = x['params']['leading_edge_radius_lower_weighting']
        max_thickness_weighting             = x['params']['max_thickness_weighting']
        max_thickness_upper_weighting       = x['params']['max_thickness_upper_weighting']
        max_thickness_lower_weighting       = x['params']['max_thickness_lower_weighting']
        radii_skew_weighting                = x['params']['radii_skew_weighting']
        curvature_weighting                 = x['params']['curvature_weighting']
        lower_surface_curvature_weighting   = x['params']['lower_surface_curvature_weighting']
        te_cone_violation_weighting         = x['params']['te_cone_violation_weighting']
        CL_target_weighting                 = x['params']['CL_target_weighting']
        clean_moment_weighting              = x['params']['clean_moment_weighting']
        rough_moment_weighting              = x['params']['rough_moment_weighting']
        cpmin_design_weighting              = x['params']['cp_min_design_weighting']
        cpmin_at_alpha_offset_weighting     = x['params']['cp_min_at_alpha_offset_weighting']
        cpmin_prestall_weighting            = x['params']['cp_min_prestall_weighting']
        infeasibility_penalty               = x['params']['infeasibility_penalty']
        min_radius_location_upper_weighting = x['params']['min_radius_location_upper_weighting']
        min_radius_location_lower_weighting = x['params']['min_radius_location_lower_weighting']
        toothpick_weighting                 = x['params']['toothpick_weighting']

        N_crit_clean     = x['params']['N_crit_clean']
        N_crit_rough     = x['params']['N_crit_rough']
        alpha_min_clean  = x['params']['alpha_min_clean']
        alpha_max_clean  = x['params']['alpha_max_clean']
        alpha_step_clean = x['params']['alpha_step_clean']
        alpha_min_rough  = x['params']['alpha_min_rough']
        alpha_max_rough  = x['params']['alpha_max_rough']
        alpha_step_rough = x['params']['alpha_step_rough']
        xtp_u_clean      = x['params']['xtp_u_clean']
        xtp_l_clean      = x['params']['xtp_l_clean']
        xtp_u_rough      = x['params']['xtp_u_rough']
        xtp_l_rough      = x['params']['xtp_l_rough']
        N_points_moi     = x['params']['N_points_moi']

        selected_tool = x['params']['tool']

        if selected_tool == 'xfoil':
            path_to_XFOIL = x['params']['xfoil_path']
            tfpre = x['params']['xfoil_tempfile_path_leader']
            xfoil_timelimit  = x['params']['xfoil_timelimit']
        elif selected_tool == 'neuralfoil':
            neuralfoil_model = x['params']['neuralfoil_model']
        else:
            raise RuntimeError('Invalid tool selection')

        # ----------------------
        # reject too large coefficients
        # ----------------------
        if max(abs(np.array(K_upper)))>2.0:
            return [pid, np.inf, np.inf, False, -20] + [0]*N_reported + [0]*N_constraints
        if max(abs(np.array(K_lower)))>2.0:
            return [pid, np.inf, np.inf, False, -30] + [0]*N_reported + [0]*N_constraints

        # ----------------------
        # build airfoil
        # ----------------------
        afl_geo = Kulfan(TE_gap = te_gap)
        afl_geo.upperCoefficients = K_upper
        afl_geo.lowerCoefficients = K_lower
        afl_geo.chord = 1.0*units.m
        
        # ----------------------
        # reject self intersecting
        # ----------------------
        hts = afl_geo.getNormalizedHeight()
        if any(hts<0):
            return [pid, np.inf, np.inf, False, -10] + [0]*N_reported + [0]*N_constraints

        # ----------------------
        # Run Clean Data
        # ----------------------
        if selected_tool == 'xfoil':
            res1 = run_xfoil('alfa', K_upper, K_lower, [alpha_min_clean, alpha_max_clean, alpha_step_clean], Re=Re, N_crit=N_crit_clean, xtp_u=xtp_u_clean, xtp_l=xtp_l_clean, TE_gap = te_gap, timelimit=xfoil_timelimit, path_to_XFOIL=path_to_XFOIL, tfpre=tfpre)
        elif selected_tool == 'neuralfoil':
            res1 = run_neuralfoil('alfa', K_upper, K_lower, [alpha_min_clean, alpha_max_clean, alpha_step_clean], Re=Re, N_crit=N_crit_clean, xtp_u=xtp_u_clean, xtp_l=xtp_l_clean, TE_gap = te_gap, model = neuralfoil_model)
        else:
            raise RuntimeError('Invalid tool selection')              

        cl_clean    = res1['cl']
        cd_clean    = res1['cd']
        cm_clean    = res1['cm']
        cpmin_clean = res1['cpmin']
        alpha_clean = res1['alpha']
        LoD_clean   = [cl_clean[i]/cd_clean[i] for i in range(0,len(cl_clean))] 

        # ----------------------
        # Run Rough Data
        # ----------------------  
        if selected_tool == 'xfoil':
            res2 = run_xfoil('alfa', K_upper, K_lower, [alpha_min_rough, alpha_max_rough, alpha_step_rough], Re=Re, N_crit=N_crit_rough, xtp_u=xtp_u_rough, xtp_l=xtp_l_rough, TE_gap = te_gap, timelimit=xfoil_timelimit, path_to_XFOIL=path_to_XFOIL, tfpre=tfpre)
        elif selected_tool == 'neuralfoil':
            res2 = run_neuralfoil('alfa', K_upper, K_lower, [alpha_min_rough, alpha_max_rough, alpha_step_rough], Re=Re, N_crit=N_crit_rough, xtp_u=xtp_u_rough, xtp_l=xtp_l_rough, TE_gap = te_gap, model = neuralfoil_model)
        else:
            raise RuntimeError('Invalid tool selection')              

        if res2 is None:
            raise RuntimeError("Xfoil failed")                

        cl_rough    = res2['cl']
        cd_rough    = res2['cd']
        cm_rough    = res2['cm']
        cpmin_rough = res2['cpmin']
        alpha_rough = res2['alpha']
        LoD_rough   = [cl_rough[i]/cd_rough[i] for i in range(0,len(cl_rough))] 
        
        # ----------------------
        # Find the stall locations
        # ----------------------
        try:
            mid_idx_clean = np.argmin(abs(np.array(alpha_clean)))
            mid_idx_rough = np.argmin(abs(np.array(alpha_rough)))

            positive_peak_index_clean = mid_idx_clean
            for i in range(mid_idx_clean, mid_idx_clean+len(alpha_clean)):
                if cl_clean[i] > cl_clean[i-1]:
                    pass
                else:
                    positive_peak_index_clean = i-1
                    break
                # will reach an index error if no peak is present

            positive_peak_index_rough = mid_idx_rough
            for i in range(mid_idx_rough, mid_idx_rough+len(alpha_rough)):
                if cl_rough[i] > cl_rough[i-1]:
                    pass
                else:
                    positive_peak_index_rough = i-1
                    break
                # will reach an index error if no peak is present

        except:
            # threw an index error, so no peak is present
            return [pid, np.inf, np.inf, False, -70] + [0]*N_reported + [0]*N_constraints
            
        # ----------------------
        # Find alpha_design
        # ----------------------

        if cl_clean[positive_peak_index_clean] > cl_design:
            alpha_design = np.interp(cl_design, 
                                     np.array(cl_clean)[0:positive_peak_index_clean], 
                                     np.array(alpha_clean)[0:positive_peak_index_clean] )
        else:
            # raise ValueError("airfoil cannot reach target CL")
            return [pid, np.inf, np.inf, False, -80] + [0]*N_reported + [0]*N_constraints
        
        # ----------------------
        # Begin objective function and constraints
        # ----------------------
        conpen = 0.0
        cons = []

        # ----------------------
        # Pure L/D performance
        # ----------------------
        LoD_clean_at_design_alpha = np.interp(alpha_design, 
                                     np.array(alpha_clean)[0:positive_peak_index_clean], 
                                     np.array(LoD_clean)[0:positive_peak_index_clean] )

        LoD_rough_at_design_alpha = np.interp(alpha_design, 
                                     np.array(alpha_rough)[0:positive_peak_index_rough], 
                                     np.array(LoD_rough)[0:positive_peak_index_rough] )
        
        obj1 = -1*LoD_clean_at_design_alpha

        obj2 = -1*LoD_rough_at_design_alpha

        # ----------------------
        # Stall margin
        # ----------------------
        stall_margin_clean = alpha_clean[positive_peak_index_clean] - alpha_design
        stall_margin_rough = alpha_rough[positive_peak_index_rough] - alpha_design
        
        # target_stall_margin = 4.0
        # 0 if valid, otherwise 4.0-target_stall_margin
        if stall_margin_clean_weighting is not None and stall_margin_clean_weighting != 0.0:
            cons.append(target_stall_margin-(min([target_stall_margin,stall_margin_clean])))
            conpen += stall_margin_clean_weighting*(target_stall_margin-min([target_stall_margin, stall_margin_clean]))
        else:
            cons.append(0.0)
        if stall_margin_rough_weighting is not None and stall_margin_rough_weighting != 0.0:
            cons.append(target_stall_margin-(min([target_stall_margin,stall_margin_rough])))
            conpen += stall_margin_rough_weighting*(target_stall_margin-min([target_stall_margin, stall_margin_rough]))
        else:
            cons.append(0.0)
        
        # ----------------------
        # lift margin
        # ----------------------
        lift_margin_clean = cl_clean[positive_peak_index_clean] - cl_design
        
        if lift_margin_clean_weighting is not None and lift_margin_clean_weighting != 0.0:
            if lift_margin_clean >=0:
                # expected behavior
                conpen += lift_margin_clean_weighting*lift_margin_clean
            else:
                # in some weird near-stall location, penalize extra
                conpen += 100*lift_margin_clean_weighting*abs(lift_margin_clean)

        
        # ----------------------
        # limit CL max
        # ----------------------
        if cl_max_limit_clean is not None and cl_max_limit_clean_weighting is not None and cl_max_limit_clean_weighting != 0.0:
            cl_max_violation_clean = cl_clean[positive_peak_index_clean] - cl_max_limit_clean
            if cl_max_violation_clean > 0:
                conpen += cl_max_limit_clean_weighting * cl_max_violation_clean
                cons.append(cl_max_violation_clean)
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)


        if cl_max_limit_rough is not None and cl_max_limit_rough_weighting is not None and cl_max_limit_rough_weighting != 0.0:
            cl_max_violation_rough = cl_rough[positive_peak_index_rough] - cl_max_limit_rough
            if cl_max_violation_rough > 0:
                conpen += cl_max_limit_rough_weighting * cl_max_violation_rough
                cons.append(cl_max_violation_rough)
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)

        # ----------------------
        # change in cl at alpha design due to roughness
        # ----------------------
        delta_cl_clean_to_rough_at_alpha_design = cl_design - np.interp(alpha_design, 
                                                                        np.array(alpha_rough)[0:positive_peak_index_rough], 
                                                                        np.array(cl_rough)[0:positive_peak_index_rough] )

        percent_delta_cl_clean_to_rough_at_alpha_design = delta_cl_clean_to_rough_at_alpha_design/cl_design
        
        # penalize if greater than 10%
        # assert(percent_delta_cl_clean_to_rough_at_alpha_design >=0)
        if delta_cl_from_roughness_weighting is not None and delta_cl_from_roughness_weighting != 0.0:
            conpen += delta_cl_from_roughness_weighting * (max([percent_delta_cl_from_roughness_threshold, abs(percent_delta_cl_clean_to_rough_at_alpha_design)])-percent_delta_cl_from_roughness_threshold)
        
        # ----------------------
        # clean L/D curve fall off
        # ----------------------
        LoD_clean_1degree_left = np.interp(alpha_design - alpha_falloff_offset, 
                                             np.array(alpha_clean)[0:positive_peak_index_clean], 
                                             np.array(LoD_clean)[0:positive_peak_index_clean] )

        LoD_clean_1degree_right = np.interp(alpha_design + alpha_falloff_offset, 
                                             np.array(alpha_clean)[0:positive_peak_index_clean], 
                                             np.array(LoD_clean)[0:positive_peak_index_clean] )

        percent_change_LoD_clean_left  = (LoD_clean_1degree_left-LoD_clean_at_design_alpha)/LoD_clean_at_design_alpha
        percent_change_LoD_clean_right = (LoD_clean_1degree_right-LoD_clean_at_design_alpha)/LoD_clean_at_design_alpha
     
        # penalize if greater than 15%
        # assert(percent_delta_LoD_clean_to_rough_at_alpha_design >=0)
        if LoD_falloff_weighting is not None and LoD_falloff_weighting != 0.0:
            conpen += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_clean_left )])-percent_LoD_falloff_threshold)
            conpen += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_clean_right)])-percent_LoD_falloff_threshold)

        # ----------------------
        # rough L/D curve fall off
        # ----------------------
        LoD_rough_1degree_left = np.interp(alpha_design - alpha_falloff_offset, 
                                             np.array(alpha_rough)[0:positive_peak_index_rough], 
                                             np.array(LoD_rough)[0:positive_peak_index_rough] )

        LoD_rough_1degree_right = np.interp(alpha_design + alpha_falloff_offset, 
                                             np.array(alpha_rough)[0:positive_peak_index_rough], 
                                             np.array(LoD_rough)[0:positive_peak_index_rough] )

        percent_change_LoD_rough_left  = (LoD_rough_1degree_left-LoD_rough_at_design_alpha)/LoD_rough_at_design_alpha
        percent_change_LoD_rough_right = (LoD_rough_1degree_right-LoD_rough_at_design_alpha)/LoD_rough_at_design_alpha

        # penalize if greater than 15%
        # assert(percent_delta_LoD_clean_to_rough_at_alpha_design >=0)
        if LoD_falloff_weighting is not None and LoD_falloff_weighting != 0.0:
            conpen += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_rough_left )])-percent_LoD_falloff_threshold)
            conpen += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_rough_right)])-percent_LoD_falloff_threshold)
       
        # ----------------------
        # structure surrogates
        # ----------------------
        moi_afl = copy.deepcopy(afl_geo)
        moi_afl.utility.Npoints = N_points_moi
        Ixx = moi_afl.Ixx.magnitude
        Iyy = moi_afl.Iyy.magnitude
        Izz = moi_afl.Izz.magnitude
        A   = moi_afl.area.magnitude
        
        Ixx_target = Ixx_con 
        Iyy_target = Iyy_con 
        Izz_target = Izz_con 
        A_target   = A_con

        # 0 if valid, otherwise target-val
        if ixx_weighting is not None and ixx_weighting != 0.0:
            cons.append(Ixx_target-(min([Ixx_target,Ixx])))
            conpen += ixx_weighting * (Ixx_target-min([Ixx_target, Ixx]))
        else:
            cons.append(0.0)
        
        if iyy_weighting is not None and iyy_weighting != 0.0:
            cons.append(Iyy_target-(min([Iyy_target,Iyy])))
            conpen += iyy_weighting * (Iyy_target-min([Iyy_target, Iyy]))
        else:
            cons.append(0.0)
        
        if izz_weighting is not None and izz_weighting != 0.0:
            cons.append(Izz_target-(min([Izz_target,Izz])))
            conpen += izz_weighting * (Izz_target-min([Izz_target, Izz]))
        else:
            cons.append(0.0)
        
        if a_weighting is not None and a_weighting != 0.0:
            cons.append(A_target-(min([A_target,A])))
            conpen += a_weighting * (A_target-min([A_target, A]))
        else:
            cons.append(0.0)
        
        # ----------------------
        # reject leading edges with too tight radii
        # ----------------------
        leading_edge_radius_upper, leading_edge_radius_lower = afl_geo.leadingEdgeRadius()
        
        leruViolation = -1*min([0,leading_edge_radius_upper-ler_con])
        lerlViolation = -1*min([0,leading_edge_radius_lower-ler_con])
        
        if leading_edge_radius_upper_weighting is not None and leading_edge_radius_upper_weighting != 0.0:
            cons.append(leruViolation)
            conpen += leading_edge_radius_upper_weighting * leruViolation
        else:
            cons.append(0.0)
        if leading_edge_radius_lower_weighting is not None and leading_edge_radius_lower_weighting != 0.0:
            cons.append(lerlViolation)
            conpen += leading_edge_radius_lower_weighting * lerlViolation
        else:
            cons.append(0.0)

        # ----------------------
        # reject violations of TE_cone
        # ----------------------
        
        height_upper_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaUpper)
        height_lower_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaLower)

        midpoint_at_98_percent = ( height_upper_at_98percent + height_lower_at_98percent )/2

        upper_10deg_cone = midpoint_at_98_percent + np.tan(cone_angle/2/180*np.pi)*(1-te_frac)
        lower_10deg_cone = midpoint_at_98_percent - np.tan(cone_angle/2/180*np.pi)*(1-te_frac)

        teViolation = 0
        for i,psi_val in enumerate(afl_geo.psi):
            if psi_val >= te_frac:
                zeta_upper_val = afl_geo.zetaUpper[i]
                zeta_lower_val = afl_geo.zetaLower[i]

                h_upper = upper_10deg_cone - (psi_val-te_frac)*(upper_10deg_cone/(1-te_frac)) + afl_geo.zetaUpper[-1]
                h_lower = lower_10deg_cone - (psi_val-te_frac)*(lower_10deg_cone/(1-te_frac)) + afl_geo.zetaLower[-1]

                if zeta_upper_val < h_upper:
                    teViolation += (h_upper-zeta_upper_val)                
                if zeta_lower_val > h_lower:
                    teViolation += (zeta_lower_val-h_lower)
        
        # Account for a small numerical issues
        if teViolation < 1e-8:
            teViolation = 0

        if te_cone_violation_weighting is not None and te_cone_violation_weighting != 0.0:
            cons.append(teViolation)
            conpen += te_cone_violation_weighting*teViolation
        else:
            cons.append(0.0)
       
        # ----------------------
        # reject if the max thickness is further forward than 25%
        # ----------------------
        tau_loc = afl_geo.taumax_psi
        if max_thickness_weighting is not None and max_thickness_weighting != 0.0 and tau_loc < max_thickness_loc:
            cons.append(max_thickness_loc-tau_loc)
            conpen += max_thickness_weighting * (max_thickness_loc-tau_loc) 
        else:
            cons.append(0)

        # ----------------------
        # reject if the max thickness on upper surface is further forward than 25%
        # ----------------------
        tau_loc_u = afl_geo.taumax_psi_upper
        if max_thickness_upper_weighting is not None and max_thickness_upper_weighting != 0.0 and tau_loc_u < max_thickness_loc_upper:
            cons.append(max_thickness_loc_upper-tau_loc_u)
            conpen += max_thickness_upper_weighting * (max_thickness_loc_upper-tau_loc_u) 
        else:
            cons.append(0)

        # ----------------------
        # reject if the max thickness on lower surface is further forward than 25%
        # ----------------------
        tau_loc_l = afl_geo.taumax_psi_lower
        if max_thickness_lower_weighting is not None and max_thickness_lower_weighting != 0.0 and tau_loc_l < max_thickness_loc_lower:
            cons.append(max_thickness_loc_lower-tau_loc_l)
            conpen += max_thickness_lower_weighting * (max_thickness_loc_lower-tau_loc_l) 
        else:
            cons.append(0)

        # ----------------------
        # attempt to stop significantly asymmetric leading edges, should not be true as this is very conservative
        # ----------------------
        radii = afl_geo.leadingEdgeRadius()
        if radii_skew_weighting is not None and radii_skew_weighting != 0.0 and max(radii) > ler_skew_factor*min(radii):
            conpen += radii_skew_weighting * (max(radii)-ler_skew_factor*min(radii))
            cons.append(max(radii)-ler_skew_factor*min(radii))
        else:
            cons.append(0)

        # ----------------------
        # Verify thickness, though this should be true by construction (except when nans were introduced)
        # ----------------------
        # negate because 0 is an unviolated constraint
        cons.append(int(not abs(afl_geo.tau-tau)<1e-4))

        # ----------------------
        # Penalize upper surface concavity (positive curvature)
        # ----------------------
        delta_zeta_upper = afl_geo.zetaUpper[1:] - afl_geo.zetaUpper[0:-1]
        delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        first_derivative_approx = (delta_zeta_upper / delta_psi)
        delta_first_derivative = first_derivative_approx[1:] - first_derivative_approx[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx = delta_first_derivative/delta_delta_psi
        positive_curvature = []
        for i in range(0, len(second_derivative_approx)):
            if second_derivative_approx[i] >0:
                positive_curvature.append(second_derivative_approx[i])
        if curvature_weighting is not None and curvature_weighting != 0.0:
            conpen += curvature_weighting * sum(positive_curvature)
            cons.append(sum(positive_curvature))
        else:
            cons.append(0.0)
        # ----------------------
        # Penalize the aft curvature
        excess_curvature = second_derivative_approx - curvature_bound
        ec_psi = afl_geo.psi[1:-1] # approximate, but close
        ec_sum = 0
        if curvature_weighting is not None and curvature_weighting != 0.0:
            for tst_idx, ecpsi in enumerate(ec_psi):
                if ecpsi >= ec_cutoff:
                    if excess_curvature[tst_idx] < 0:
                        ec_sum += curvature_weighting * abs(excess_curvature[tst_idx])
            conpen += ec_sum
            cons.append(ec_sum/curvature_weighting)
        else:
            cons.append(0.0)

        # ----------------------
        # Penalize lower surface curvature changes if it happens more than once
        # ----------------------
        delta_zeta_lower = afl_geo.zetaLower[1:] - afl_geo.zetaLower[0:-1]
        delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        first_derivative_approx_l = (delta_zeta_lower / delta_psi)
        delta_first_derivative_l = first_derivative_approx_l[1:] - first_derivative_approx_l[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx_l = delta_first_derivative_l/delta_delta_psi
        sflips = 0
        sgn = second_derivative_approx_l[0]/abs(second_derivative_approx_l[0])
        for i in range(0, len(second_derivative_approx_l)):
            if second_derivative_approx_l[i]/abs(second_derivative_approx_l[i]) != sgn:
                sgn = second_derivative_approx_l[i]/abs(second_derivative_approx_l[i])
                sflips += 1
        if lower_surface_curvature_weighting is not None and lower_surface_curvature_weighting != 0.0 and sflips > 1:
            conpen += lower_surface_curvature_weighting * (sflips-1)
            cons.append(sflips-1)
        else:
            cons.append(0)

        # ----------------------
        # Ensure 36% airfoil has sufficient CL for rough case
        # ----------------------
        if target_alpha is not None and CL_target_weighting is not None and CL_target_weighting != 0.0:
            if max(alpha_rough)<target_alpha:
                #  Higher alphas did not converge
                return [pid, np.inf, np.inf, False, -60] + [0]*N_reported + [0]*N_constraints

            cl_at_target_degrees_rough = np.interp(target_alpha, 
                                               np.array(alpha_rough)[0:positive_peak_index_rough], 
                                               np.array(cl_rough)[0:positive_peak_index_rough] )

            if cl_at_target_degrees_rough <= target_cl:
                conpen += CL_target_weighting * abs(target_cl-cl_at_target_degrees_rough)
                cn = target_cl-cl_at_target_degrees_rough
            else:
                cn = 0.0

            cons.append(cn)
        else:
            cons.append(0.0)

        # ----------------------
        # Add a moment constraint (clean)
        # ----------------------
        if CMc_in is not None and clean_moment_weighting is not None and clean_moment_weighting != 0.0:
            try:
                cm_design_clean = -1*abs(CMc_in)
                alpha_design_range_indicies = [i for i,alpha_val in enumerate(alpha_clean) if alpha_val >= alpha_design-cm_alpha_band and alpha_val <= alpha_design+cm_alpha_band]
                min_cm_clean = min([cm_clean[i] for i in alpha_design_range_indicies])
                if min_cm_clean < cm_design_clean:
                    conpen += clean_moment_weighting * abs(cm_design_clean - min_cm_clean)
                    cons.append(abs(cm_design_clean - min_cm_clean))
                else:
                    cons.append(0.0)
            except:
                # something strange happened
                conpen += clean_moment_weighting * 100
                cons.append(clean_moment_weighting)
        else:
            cons.append(0.0)
            
        # ----------------------
        # Add a moment constraint (rough)
        # ----------------------
        if CMr_in is not None and rough_moment_weighting is not None and rough_moment_weighting != 0.0:
            try:
                cm_design_rough = -1*abs(CMr_in)
                alpha_design_range_indicies = [i for i,alpha_val in enumerate(alpha_rough) if alpha_val >= alpha_design-cm_alpha_band and alpha_val <= alpha_design+cm_alpha_band]
                min_cm_rough = min([cm_rough[i] for i in alpha_design_range_indicies])
                if min_cm_rough < cm_design_rough:
                    conpen += rough_moment_weighting * abs(cm_design_rough - min_cm_rough)
                    cons.append(abs(cm_design_rough - min_cm_rough))
                else:
                    cons.append(0.0)
            except:
                # something strange happened
                conpen += rough_moment_weighting * 100
                cons.append(rough_moment_weighting)
        else:
            cons.append(0.0)

        # ----------------------
        # Add a minimum pressure coefficient constraint at design condition
        # ----------------------
        cpmin_swept_clean_design = np.interp(alpha_design, alpha_clean[0:positive_peak_index_clean], cpmin_clean[0:positive_peak_index_clean])
        cpmin_swept_rough_design = np.interp(alpha_design, alpha_rough[0:positive_peak_index_rough], cpmin_rough[0:positive_peak_index_rough])

        if cp_min_design is not None and cpmin_design_weighting is not None and cpmin_design_weighting != 0.0:
            if cpmin_swept_clean_design < cp_min_design:
                conpen += cpmin_design_weighting * abs(cp_min_design - cpmin_swept_clean_design)
                cons.append(abs(cp_min_design - cpmin_swept_clean_design))
            else:
                cons.append(0.0)

            if cpmin_swept_rough_design < cp_min_design:
                conpen += cpmin_design_weighting * abs(cp_min_design - cpmin_swept_rough_design)
                cons.append(abs(cp_min_design - cpmin_swept_rough_design))
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)
            cons.append(0.0)

        # ----------------------
        # Add a minimum pressure coefficient constraint at an offset from design condition
        # ----------------------
        
        if cp_min_at_alpha_offset is not None and cpmin_at_alpha_offset_weighting is not None and cpmin_at_alpha_offset_weighting != 0.0:
            cpmin_swept_clean_design = np.interp(alpha_design+cp_min_alpha_offset, alpha_clean[0:positive_peak_index_clean], cpmin_clean[0:positive_peak_index_clean])
            cpmin_swept_rough_design = np.interp(alpha_design+cp_min_alpha_offset, alpha_rough[0:positive_peak_index_rough], cpmin_rough[0:positive_peak_index_rough])

            if cpmin_swept_clean_design < cp_min_at_alpha_offset:
                conpen += cpmin_at_alpha_offset_weighting * abs(cp_min_at_alpha_offset - cpmin_swept_clean_design)
                cons.append(abs(cp_min_at_alpha_offset - cpmin_swept_clean_design))
            else:
                cons.append(0.0)

            if cpmin_swept_rough_design < cp_min_at_alpha_offset:
                conpen += cpmin_at_alpha_offset_weighting * abs(cp_min_at_alpha_offset - cpmin_swept_rough_design)
                cons.append(abs(cp_min_at_alpha_offset - cpmin_swept_rough_design))
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)
            cons.append(0.0)

        # ----------------------
        # Add a minimum pressure coefficient constraint pre-stall
        # ----------------------
        cpmin_swept_clean_prestall = min( np.array(cpmin_clean)[0:positive_peak_index_clean] )
        cpmin_swept_rough_prestall = min( np.array(cpmin_rough)[0:positive_peak_index_rough] )

        if cp_min_prestall is not None and cpmin_prestall_weighting is not None and cpmin_prestall_weighting != 0.0:
            if cpmin_swept_clean_prestall < cp_min_prestall:
                conpen += cpmin_prestall_weighting * abs(cp_min_prestall - cpmin_swept_clean_prestall)
                cons.append(abs(cp_min_prestall - cpmin_swept_clean_prestall))
            else:
                cons.append(0.0)

            if cpmin_swept_rough_prestall < cp_min_prestall:
                conpen += cpmin_prestall_weighting * abs(cp_min_prestall - cpmin_swept_rough_prestall)
                cons.append(abs(cp_min_prestall - cpmin_swept_rough_prestall))
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)
            cons.append(0.0)

        # ----------------------
        # Constrain the locations of the min radius of curvature
        # ----------------------
        delta_zeta_upper = afl_geo.zetaUpper[1:] - afl_geo.zetaUpper[0:-1]
        delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        first_derivative_approx = (delta_zeta_upper / delta_psi)
        first_derivative_average = (first_derivative_approx[1:] + first_derivative_approx[0:-1]) / 2.0
        delta_first_derivative = first_derivative_approx[1:] - first_derivative_approx[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx = delta_first_derivative/delta_delta_psi
        radius_of_curvature_approx =  (1+first_derivative_average**2)**1.5 / abs(second_derivative_approx)
        chopped_roc = radius_of_curvature_approx[afl_geo.psi[1:-1] <= min_radius_location_cutoff]
        computed_min_radius_loc_upper = afl_geo.psi[1:-1][afl_geo.psi[1:-1] <= min_radius_location_cutoff][np.argmin(chopped_roc)]
        if min_radius_location_upper_weighting is not None and min_radius_location_upper_weighting != 0.0 and computed_min_radius_loc_upper > min_radius_location_upper:
            conpen += min_radius_location_upper_weighting * (computed_min_radius_loc_upper - min_radius_location_upper)
            cons.append(computed_min_radius_loc_upper - min_radius_location_upper)
        else:
            cons.append(0.0)

        delta_zeta_lower = afl_geo.zetaLower[1:] - afl_geo.zetaLower[0:-1]
        delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        first_derivative_approx = (delta_zeta_lower / delta_psi)
        first_derivative_average = (first_derivative_approx[1:] + first_derivative_approx[0:-1]) / 2.0
        delta_first_derivative = first_derivative_approx[1:] - first_derivative_approx[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx = delta_first_derivative/delta_delta_psi
        radius_of_curvature_approx =  (1+first_derivative_average**2)**1.5 / abs(second_derivative_approx)
        chopped_roc = radius_of_curvature_approx[afl_geo.psi[1:-1] <= min_radius_location_cutoff]
        computed_min_radius_loc_lower = afl_geo.psi[1:-1][afl_geo.psi[1:-1] <= min_radius_location_cutoff][np.argmin(chopped_roc)]
        if min_radius_location_lower_weighting is not None and min_radius_location_lower_weighting != 0.0 and computed_min_radius_loc_lower > min_radius_location_lower:
            conpen += min_radius_location_lower_weighting * (computed_min_radius_loc_lower - min_radius_location_lower)
            cons.append(computed_min_radius_loc_lower - min_radius_location_lower)
        else:
            cons.append(0.0)


        # ----------------------
        # Toothpick constraint
        # ----------------------
        if toothpick_weighting is not None and toothpick_weighting != 0.0:
            toothpick_height_at_location = afl_geo.getNormalizedHeight(toothpick_location)
            if toothpick_height_at_location < toothpick_height:
                conpen += toothpick_weighting * (toothpick_height - toothpick_height_at_location)
                cons.append(toothpick_height - toothpick_height_at_location)
            else:
                cons.append(0.0)
        else:
            cons.append(0.0)

        # ----------------------
        # return
        # ----------------------
        con_tag = all([c==0 for c in cons])

        # Add penalty for violating the con tag
        # Helps at the end of the run
        if infeasibility_penalty is not None and infeasibility_penalty != 0.0:
            conpen += (1-con_tag) * infeasibility_penalty
        
        r_list = [
                pid,
                obj1 + conpen,
                obj2 + conpen,
                con_tag,
                alpha_design,
                LoD_clean_at_design_alpha,
                LoD_rough_at_design_alpha,
                stall_margin_clean,
                stall_margin_rough,
                lift_margin_clean,
                delta_cl_clean_to_rough_at_alpha_design,
                LoD_clean_1degree_left,
                LoD_clean_1degree_right,
                afl_geo.tau.magnitude,
                leading_edge_radius_upper, 
                leading_edge_radius_lower,
                Ixx,
                Iyy,
                Izz,
                A,
                min([cpmin_swept_clean_design, cpmin_swept_rough_design]),
            ] + cons 

        for rix, rtv in enumerate(r_list):
            if isinstance(rtv, units.Quantity):
                assert(rtv.units == units.dimensionless)
                r_list[rix] = rtv.to('dimensionless').magnitude

        return r_list
    except:
        return [pid, np.inf, np.inf, False, -90] + [0]*N_reported + [0]*N_constraints



if __name__ == '__main__':
    import time
    t0 = time.time()
    import matplotlib.pyplot as plt
    import yaml
    airfoil_seeds = [
        # [3.040459184728062092e-01, 3.865158898232732843e-01, 1.316538699269504120e-01 ,3.139820595993685348e-01 ,-2.593332834032758827e-01 ,-1.465580816167938449e-01 ,-4.868709152015524566e-01, 3.030523043601500155e-01],
        # [2.436455855485705202e-01, 5.627052809097938812e-01, 1.841844605039789085e-01 ,3.441977319654409007e-01 ,-2.291428721230603649e-01 ,-1.068366555679247792e-01 ,-2.774852502080942251e-01, 3.132586099910192323e-01],
        [.1,.1,.1,-.1,-.1,-.1,-.1,0.2],
        # [0.2895529384519609, 0.5707988120123219, 0.34602093353667224, 0.7675267435685507, 0.47442973419723866, 0.8705928599649991, 0.4336960256541951, 0.4663320196324578, -0.15104922613749502, -0.05073712004691763, -0.25913307521796636, 0.2605533024015594, -0.5160313871452309, 0.44225025302983123, 0.14956381466541888, 0.3174515752976585],
        # [2.418060954461536127e-01, 3.178580125622825769e-01, 4.933624366764263192e-01, 2.544814149655084679e-01, 6.142654843995031255e-01, 6.123154352426792846e-01, 6.045152386153840318e-01, 5.694143537925553389e-01, -1.974544732709675832e-01, -1.457661744725805009e-01, -2.960661180175324050e-06, -4.416861340205950892e-01, 3.961037349798383206e-02, 2.442436568921833862e-01, 3.393085532873446053e-01, 3.997600771488830085e-01],
        # [0.2759924824591029, 0.36505371575699025, 0.3807132511161153, 0.5050076534019967, 0.4032220194306991, 0.2769723823719758, 0.3983287124814639, 0.35820334426455686, -0.22510024324231004, -0.07854031784219763, -0.36994764134972546, -0.26914104884527795, -0.07071138196387987, 0.11279515592192352, 0.25641832721484725, 0.25250387348994047],
        # [0.27200518835149395, 0.5349785727756038, 0.08526832357342562, 0.8798051137895084, -0.013288258418964789, 0.6398114717677857, 0.19013282950501767, 0.40686735458401174, -0.19584565656757194, -0.249756231310527, -0.19629487910041066, -0.2311328308373759, -0.17125182017854618, -0.06893729120760746, 0.17150653015921413, 0.20454160490141127],
        # 
        # 
        [0.27200518835149395, 0.5349785727756038, 0.08526832357342562, 0.8798051137895084, -0.013288258418964789, 0.6398114717677857, 0.19013282950501767, 0.40686735458401174, -0.25, -0.249756231310527, -0.19629487910041066, -0.2311328308373759, -0.17125182017854618, -0.06893729120760746, 0.17150653015921413, 0.20454160490141127]
    ]

    labels = [
        'pid',
        'obj1',
        'obj2',
        'con_tag',
        'alpha_design',
        'LoD_clean_at_design',
        'LoD_rough_at_design',
        'stall_margin_clean',
        'stall_margin_rough',
        'lift_margin_clean',
        'delta_cl_from_roughness',
        'LoD_c_1d_left',
        'LoD_c_1d_right',
        'tau',
        'ler_upper' ,
        'ler_lower',
        'Ixx',
        'Iyy',
        'Izz',
        'A',
        'cpmin',
        'con_sm_clean',
        'con_sm_rough',
        'con_clmax_clean',
        'con_clmax_rough',
        'con_ixx',
        'con_iyy',
        'con_izz',
        'con_a',
        'con_leru',
        'con_lerl',
        'con_te_cone',
        'con_max_tau',
        'con_max_tau_u',
        'con_max_tau_l',
        'con_ler_skew',
        'con_tau',
        'con_concave',
        'con_aftcurve',
        'con_lower_flips',
        'con_10deg',
        'con_mom_c',
        'con_mom_r',
        'con_cpmin_design_clean',
        'con_cpmin_design_rough',
        'con_cpmin_offset_clean',
        'con_cpmin_offset_rough',
        'con_cpmin_prestall_clean',
        'con_cpmin_prestall_rough',
        'con_min_rad_loc_upper',
        'con_min_rad_loc_lower',
        'con_toothpick',
    ]

    # params = yaml.safe_load(open("example_hydro.yaml"))
    params = yaml.safe_load(open("run_hydro_003.json"))

    for i,ind in enumerate(airfoil_seeds):
        K_upper = ind[0:int(len(ind)/2)]
        K_lower = ind[int(len(ind)/2):]

        afl_geo = Kulfan(TE_gap = params['TE_gap'])
        afl_geo.upperCoefficients = K_upper
        afl_geo.lowerCoefficients = K_lower
        afl_geo.scaleThickness(params['tau'])
        ind = afl_geo.upperCoefficients.tolist() + afl_geo.lowerCoefficients.tolist()

        x = {}
        x['pid'] = i
        x['individual'] = ind
        x['params'] = params

        res = airfoil_fitness(x)
        # print(res)
        print(len(res))
        print(len(labels))
        assert(len(res)==len(labels))
        # print(time.time()-t0)
        for ii,r in enumerate(res):
            print(labels[ii].ljust(30),': ',r)
        print()


