import numpy as np
from kulfan import Kulfan, units
import copy

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()
from xfoil_wrapper_noprint import run as run_xfoil

def airfoil_fitness(x):
    # ----------------------
    # unpack
    # ----------------------
    pid = x['pid']
    K_upper = x['individual'][0:int(len(x['individual'])/2)]
    K_lower = x['individual'][int(len(x['individual'])/2):]
    tau = x['tau']

    try:
        design_matrix = [
            # tau,  CL,  spn,     Re
            [0.15, 1.5, 1.00, 10.0e6, ],
            [0.18, 1.5, 1.00, 10.0e6, ],
            [0.21, 1.5, 1.00, 12.0e6, ],
            [0.24, 1.4, 0.85, 13.0e6, ],
            [0.27, 1.3, 0.55, 16.0e6, ],
            [0.30, 1.2, 0.50, 18.0e6, ],
            [0.33, 1.2, 0.35, 16.0e6, ],
            [0.36, 1.2, 0.20, 13.0e6, ],
        ]

        te_gap_lookup = {
            '15':  0.00196,
            '18':  0.00230,
            '21':  0.00262,
            '24':  0.00751,
            '27':  0.01012,
            '30':  0.01140,
            '33':  0.01140,
            '36':  0.01140,
        }
        te_gap = te_gap_lookup[str(int(100*tau))]

        objective_terms_lookup = {
            #       cln, rt,  rm
            '15':  [230, 119, 124],
            '18':  [220, 118, 123],
            '21':  [190, 117, 122],
            '24':  [160, 113, 118],
            '27':  [155, 106, 111],
            '30':  [140,  93,  98],
            '33':  [130,  85,  90],
            '36':  [125,  73,  78],
        }

        cone_angle_lookup = {
            '15': 10.0,
            '18': 10.0,
            '21': 5.0,
            '24': 5.0,
            '27': 5.0,
            '30': 0.0,
            '33': 0.0,
            '36': 0.0,
        }

        Ixx_con = [0.00011000, 0.00017438, 0.00027518, 0.00041096, 0.00058321, 0.00079640, 0.00105795, 0.00137822 ]
        Iyy_con = [0.00397999, 0.00436351, 0.00493714, 0.00561409, 0.00633417, 0.00706380, 0.00779600, 0.00855043 ]
        Izz_con = [0.00408809, 0.00454606, 0.00521632, 0.00602287, 0.00691323, 0.00785849, 0.00885328, 0.00991577 ]
        A_con   = [0.08700496, 0.09995900, 0.11477620, 0.13051205, 0.14660942, 0.16289864, 0.17959744, 0.19731100 ]

        ler_con = [0.007, 0.008, 0.01, 0.025, 0.03, 0.04, 0.06, 0.08]

        N_reported    = 0
        N_constraints = 0
        N_points_moi  = 20

        constraint_base_penalty = 1e2

        target_stall_margin                      = 4.0
        pecent_delta_cl_from_roughness_threshold = 0.10
        alpha_falloff_offset                     = 1.0
        percent_LoD_falloff_threshold            = 0.15
        max_thickness_loc                        = 0.275
        max_thickness_loc_upper                  = 0.275
        max_thickness_loc_lower                  = 0.275
        ler_skew_factor                          = 1.9
        cone_angle                               = cone_angle_lookup[str(int(100*tau))]
        te_frac                                  = 0.95
        LoD_reward_threshold_clean               = objective_terms_lookup[str(int(100*tau))][0]
        LoD_reward_threshold_rough               = objective_terms_lookup[str(int(100*tau))][1]
        LoD_reward_threshold_max_rough           = objective_terms_lookup[str(int(100*tau))][2]

        rough_to_clean_mulFactor            = 1
        clean_LoD_weighting                 = 1
        rough_LoD_weighting                 = clean_LoD_weighting * rough_to_clean_mulFactor
        stall_margin_clean_weighting        = constraint_base_penalty
        stall_margin_rough_weighting        = constraint_base_penalty
        lift_margin_clean_weighting         = 0.5
        delta_cl_from_roughness_weighting   = 10 * constraint_base_penalty
        LoD_falloff_weighting               = 0.5 * constraint_base_penalty
        ixx_weighting                       = 1e6
        iyy_weighting                       = 1e4
        izz_weighting                       = 1e4
        a_weighting                         = 1e4
        leading_edge_radius_upper_weighting = 10*constraint_base_penalty
        leading_edge_radius_lower_weighting = 10*constraint_base_penalty
        max_thickness_weighting             = 100*constraint_base_penalty
        max_thickness_upper_weighting       = 100*constraint_base_penalty
        max_thickness_lower_weighting       = 500*constraint_base_penalty
        radii_skew_weighting                = 10*constraint_base_penalty
        curvature_weighting                 = 1
        lower_surface_curvature_weighting   = constraint_base_penalty
        te_cone_violation_weighting         = 1000*constraint_base_penalty
        flareout_weighting                  = 100*constraint_base_penalty
        
        rough_constraint_penalty  = 10*constraint_base_penalty
        # ----------------------
        # reject too large coefficients
        # ----------------------
        if max(abs(np.array(K_upper)))>2.0:
            return [pid, np.inf, False, -20] + [0]*N_reported + [0]*N_constraints
        if max(abs(np.array(K_lower)))>2.0:
            return [pid, np.inf, False, -30] + [0]*N_reported + [0]*N_constraints

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
            return [pid, np.inf, False, -10] + [0]*N_reported + [0]*N_constraints

        # ----------------------
        # Run Clean Data
        # ----------------------
        cl_design = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][1]
        Re = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][3]
        
        fastrun_clean_success = False
        res1 = run_xfoil('alfa', K_upper, K_lower, [0, 30,1], Re=Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap, timelimit=15)
        if res1 is None:
            raise RuntimeError("Xfoil failed")                

        cl_clean    = res1['cl']
        cd_clean    = res1['cd']
        alpha_clean = res1['alpha']
        LoD_clean   = [cl_clean[i]/cd_clean[i] for i in range(0,len(cl_clean))] 

        # ----------------------
        # Run Rough Data
        # ----------------------  
        res2 = run_xfoil('alfa', K_upper, K_lower, [0, 20,1], Re=Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap, timelimit=15)
        if res2 is None:
            raise RuntimeError("Xfoil failed")                

        cl_rough    = res2['cl']
        cd_rough    = res2['cd']
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
            return [pid, np.inf, False, -70] + [0]*N_reported + [0]*N_constraints
            
        # ----------------------
        # Find alpha_design
        # ----------------------
        cl_design = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][1]

        if cl_clean[positive_peak_index_clean] > cl_design:
            alpha_design = np.interp(cl_design, 
                                     np.array(cl_clean)[0:positive_peak_index_clean], 
                                     np.array(alpha_clean)[0:positive_peak_index_clean] )
        else:
            # raise ValueError("airfoil cannot reach target CL")
            return [pid, np.inf, False, -80] + [0]*N_reported + [0]*N_constraints
        
        # ----------------------
        # Begin objective function and constraints
        # ----------------------
        obj  = 0.0
        cons = []
        con_tag = True

        # ----------------------
        # Pure L/D performance
        # ----------------------
        LoD_clean_at_design_alpha = np.interp(alpha_design, 
                                     np.array(alpha_clean)[0:positive_peak_index_clean], 
                                     np.array(LoD_clean)[0:positive_peak_index_clean] )

        LoD_rough_at_design_alpha = np.interp(alpha_design, 
                                     np.array(alpha_rough)[0:positive_peak_index_rough], 
                                     np.array(LoD_rough)[0:positive_peak_index_rough] )
        
        obj += clean_LoD_weighting*-1*(LoD_clean_at_design_alpha-LoD_reward_threshold_clean)

        if LoD_rough_at_design_alpha > LoD_reward_threshold_max_rough :
            obj += rough_LoD_weighting*-1*(LoD_rough_at_design_alpha-LoD_reward_threshold_max_rough)
        else:
            obj += rough_constraint_penalty*rough_LoD_weighting*-1*(LoD_rough_at_design_alpha-LoD_reward_threshold_max_rough)
        
        # ----------------------
        # Stall margin
        # ----------------------
        stall_margin_clean = alpha_clean[positive_peak_index_clean] - alpha_design
        stall_margin_rough = alpha_rough[positive_peak_index_rough] - alpha_design
        
        # target_stall_margin = 4.0
        # 0 if valid, otherwise 4.0-target_stall_margin
        cons.append(target_stall_margin-(min([target_stall_margin,stall_margin_clean])))
        cons.append(target_stall_margin-(min([target_stall_margin,stall_margin_rough])))
        
        obj += stall_margin_clean_weighting*(target_stall_margin-min([target_stall_margin, stall_margin_clean])) 
        obj += stall_margin_rough_weighting*(target_stall_margin-min([target_stall_margin, stall_margin_rough])) 
        
        # ----------------------
        # lift margin
        # ----------------------
        lift_margin_clean = cl_clean[positive_peak_index_clean] - cl_design
        
        if lift_margin_clean >=0:
            # expected behavior
            obj += lift_margin_clean_weighting*lift_margin_clean
        else:
            # in some weird near-stall location, penalize
            obj += lift_margin_clean_weighting*abs(lift_margin_clean)*constraint_base_penalty

        # ----------------------
        # change in cl at alpha design due to roughness
        # ----------------------
        delta_cl_clean_to_rough_at_alpha_design = cl_design - np.interp(alpha_design, 
                                                                        np.array(alpha_rough)[0:positive_peak_index_rough], 
                                                                        np.array(cl_rough)[0:positive_peak_index_rough] )

        percent_delta_cl_clean_to_rough_at_alpha_design = delta_cl_clean_to_rough_at_alpha_design/cl_design
        
        # penalize if greater than 10%
        # assert(percent_delta_cl_clean_to_rough_at_alpha_design >=0)
        obj += delta_cl_from_roughness_weighting * (max([pecent_delta_cl_from_roughness_threshold, abs(percent_delta_cl_clean_to_rough_at_alpha_design)])-pecent_delta_cl_from_roughness_threshold)
        
        # ----------------------
        # clean L/D curve fall off
        # ----------------------
        # LoD_clean_1degree_left = np.interp(0.85*alpha_design, 
        LoD_clean_1degree_left = np.interp(alpha_design - alpha_falloff_offset, 
                                             np.array(alpha_clean)[0:positive_peak_index_clean], 
                                             np.array(LoD_clean)[0:positive_peak_index_clean] )

        # LoD_clean_1degree_right = np.interp(1.15*alpha_design, 
        LoD_clean_1degree_right = np.interp(alpha_design + alpha_falloff_offset, 
                                             np.array(alpha_clean)[0:positive_peak_index_clean], 
                                             np.array(LoD_clean)[0:positive_peak_index_clean] )

        percent_change_LoD_clean_left  = (LoD_clean_1degree_left-LoD_clean_at_design_alpha)/LoD_clean_at_design_alpha
        percent_change_LoD_clean_right = (LoD_clean_1degree_right-LoD_clean_at_design_alpha)/LoD_clean_at_design_alpha
     
        # penalize if greater than 15%
        # assert(percent_delta_LoD_clean_to_rough_at_alpha_design >=0)
        obj += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_clean_left )])-percent_LoD_falloff_threshold)
        obj += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_clean_right)])-percent_LoD_falloff_threshold)

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
        obj += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_rough_left )])-percent_LoD_falloff_threshold)
        obj += LoD_falloff_weighting * (max([percent_LoD_falloff_threshold, abs(percent_change_LoD_rough_right)])-percent_LoD_falloff_threshold)
       
        # ----------------------
        # structure surrogates
        # ----------------------
        moi_afl = copy.deepcopy(afl_geo)
        moi_afl.utility.Npoints = N_points_moi  = 20
        Ixx = moi_afl.Ixx.magnitude
        Iyy = moi_afl.Iyy.magnitude
        Izz = moi_afl.Izz.magnitude
        A   = moi_afl.area.magnitude
        
        Ixx_target = Ixx_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        Iyy_target = Iyy_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        Izz_target = Izz_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        A_target   = A_con[[dmr[0] for dmr in design_matrix].index(tau)] 

        # 0 if valid, otherwise target-val
        cons.append(Ixx_target-(min([Ixx_target,Ixx])))
        obj += ixx_weighting * (Ixx_target-min([Ixx_target, Ixx]))
        
        cons.append(Iyy_target-(min([Iyy_target,Iyy])))
        obj += iyy_weighting * (Iyy_target-min([Iyy_target, Iyy]))
        
        cons.append(Izz_target-(min([Izz_target,Izz])))
        obj += izz_weighting * (Izz_target-min([Izz_target, Izz])) 
        
        cons.append(A_target-(min([A_target,A])))
        obj += a_weighting * (A_target-min([A_target, A])) 
        
        # ----------------------
        # reject leading edges with too tight radii
        # ----------------------
        leading_edge_radius_upper, leading_edge_radius_lower = afl_geo.leadingEdgeRadius()
        
        leruViolation = -1*min([0,leading_edge_radius_upper-ler_con[[dmr[0] for dmr in design_matrix].index(tau)]])
        lerlViolation = -1*min([0,leading_edge_radius_lower-ler_con[[dmr[0] for dmr in design_matrix].index(tau)]])
        
        cons.append(leruViolation)
        cons.append(lerlViolation)
        
        obj += leading_edge_radius_upper_weighting * leruViolation 
        obj += leading_edge_radius_lower_weighting * lerlViolation

        # # ----------------------
        # # reject violations of TE_cone
        # # ----------------------
        # height_upper_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaUpper)
        # height_lower_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaLower)

        # midpoint_at_98_percent = ( height_upper_at_98percent + height_lower_at_98percent )/2

        # upper_10deg_cone = midpoint_at_98_percent + np.tan(cone_angle/2/180*np.pi)*(1-te_frac)
        # lower_10deg_cone = midpoint_at_98_percent - np.tan(cone_angle/2/180*np.pi)*(1-te_frac)

        # teViolation = 0
        # for i,psi_val in enumerate(afl_geo.psi):
        #     if psi_val >= te_frac:
        #         zeta_upper_val = afl_geo.zetaUpper[i]
        #         zeta_lower_val = afl_geo.zetaLower[i]

        #         h_upper = upper_10deg_cone - (psi_val-te_frac)*(upper_10deg_cone/(1-te_frac)) + afl_geo.zetaUpper[-1]
        #         h_lower = lower_10deg_cone - (psi_val-te_frac)*(lower_10deg_cone/(1-te_frac)) + afl_geo.zetaLower[-1]

        #         if zeta_upper_val < h_upper:
        #             teViolation += (h_upper-zeta_upper_val)                
        #         if zeta_lower_val > h_lower:
        #             teViolation += (zeta_lower_val-h_lower)
        
        # # Account for a small numerical issues
        # if teViolation < 1e-8:
        #     teViolation = 0

        # cons.append(teViolation)
        # obj += te_cone_violation_weighting*teViolation         
       
        # # ----------------------
        # # reject if the max thickness is further forward than 25%
        # # ----------------------
        # tau_loc = afl_geo.taumax_psi
        # if tau_loc < max_thickness_loc:
        #     cons.append(max_thickness_loc-tau_loc)
        #     obj += max_thickness_weighting * (max_thickness_loc-tau_loc) 
        # else:
        #     cons.append(0)

        # # ----------------------
        # # reject if the max thickness on upper surface is further forward than 25%
        # # ----------------------
        # tau_loc_u = afl_geo.taumax_psi_upper
        # if tau_loc_u < max_thickness_loc_upper:
        #     cons.append(max_thickness_loc_upper-tau_loc_u)
        #     obj += max_thickness_upper_weighting * (max_thickness_loc_upper-tau_loc_u) 
        # else:
        #     cons.append(0)

        # # ----------------------
        # # reject if the max thickness on upper surface is further forward than 25%
        # # ----------------------
        # tau_loc_l = afl_geo.taumax_psi_lower
        # if tau_loc_l < max_thickness_loc_lower:
        #     cons.append(max_thickness_loc_lower-tau_loc_l)
        #     obj += max_thickness_lower_weighting * (max_thickness_loc_lower-tau_loc_l) 
        # else:
        #     cons.append(0)

        # # ----------------------
        # # attempt to stop significantly asymmetric leading edges, should not be true as this is very conservative
        # # ----------------------
        # radii = afl_geo.leadingEdgeRadius()
        # if max(radii) > ler_skew_factor*min(radii):
        #     obj += radii_skew_weighting * (max(radii)-ler_skew_factor*min(radii))
        #     cons.append(max(radii)-ler_skew_factor*min(radii))
        # else:
        #     cons.append(0)

        # # ----------------------
        # # Verify thickness, though this should be true by construction (except when nans were introduced)
        # # ----------------------
        # # negate because 0 is an unviolated constraint
        # cons.append(int(not abs(afl_geo.tau-tau)<1e-4))

        # # ----------------------
        # # Penalize upper surface concavity (positive curvature)
        # # ----------------------
        # delta_zeta_upper = afl_geo.zetaUpper[1:] - afl_geo.zetaUpper[0:-1]
        # delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        # first_derivative_approx = (delta_zeta_upper / delta_psi)
        # delta_first_derivative = first_derivative_approx[1:] - first_derivative_approx[0:-1]
        # delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        # second_derivative_approx = delta_first_derivative/delta_delta_psi
        # positive_curvature = []
        # for i in range(0, len(second_derivative_approx)):
        #     if second_derivative_approx[i] >0:
        #         positive_curvature.append(second_derivative_approx[i])
        # obj += curvature_weighting * sum(positive_curvature)
        # cons.append(sum(positive_curvature))

        # # ----------------------
        # # Penalize lower surface curvature changes if it happens more than once
        # # ----------------------
        # delta_zeta_lower = afl_geo.zetaLower[1:] - afl_geo.zetaLower[0:-1]
        # delta_psi = afl_geo.psi[1:] - afl_geo.psi[0:-1]
        # first_derivative_approx_l = (delta_zeta_lower / delta_psi)
        # delta_first_derivative_l = first_derivative_approx_l[1:] - first_derivative_approx_l[0:-1]
        # delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        # second_derivative_approx_l = delta_first_derivative_l/delta_delta_psi
        # sflips = 0
        # sgn = second_derivative_approx_l[0]/abs(second_derivative_approx_l[0])
        # for i in range(0, len(second_derivative_approx_l)):
        #     if second_derivative_approx_l[i]/abs(second_derivative_approx_l[i]) != sgn:
        #         sgn = second_derivative_approx_l[i]/abs(second_derivative_approx_l[i])
        #         sflips += 1
        # if sflips > 1:
        #     obj += lower_surface_curvature_weighting * (sflips-1)
        #     cons.append(sflips-1)
        # else:
        #     cons.append(0)

        # # ----------------------
        # # Prevent a flat back style flare out
        # # ----------------------
        # nmhts = afl_geo.getNormalizedHeight()
        # maxix = np.argmax(nmhts)
        # if any(nmhts[maxix:]<nmhts[-1]):
        #     obj += (nmhts[-1]-min(nmhts[maxix:-1]))*flareout_weighting
        #     cons.append(nmhts[-1]-min(nmhts[maxix:-1]))
        # else:
        #     cons.append(0)

        # ----------------------
        # return
        # ----------------------
        con_tag = all([c==0 for c in cons])
        
        return [pid, obj, con_tag, 0.0]

        # return [
        #         pid,
        #         obj,
        #         con_tag,
        #         alpha_design,
        #         LoD_clean_at_design_alpha,
        #         LoD_rough_at_design_alpha,
        #         stall_margin_clean,
        #         stall_margin_rough,
        #         lift_margin_clean,
        #         delta_cl_clean_to_rough_at_alpha_design,
        #         LoD_clean_1degree_left,
        #         LoD_clean_1degree_right,
        #         afl_geo.tau.magnitude,
        #         leading_edge_radius_upper, 
        #         leading_edge_radius_lower,
        #         Ixx,
        #         Iyy,
        #         Izz,
        #         A,
        #     ] + cons 
    except:
        return [pid, np.inf, False, -90] + [0]*N_reported + [0]*N_constraints

