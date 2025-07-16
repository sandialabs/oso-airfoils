import numpy as np
from kulfan import Kulfan, units
from xfoil_wrapper import run as run_xfoil



def tryFastrun(K_upper, K_lower, cl_design, Re, N_crit, xtr_u, xtr_l, TE_gap):
    # try fastrun
    fastrunWorks = False

    res_fastrun = []

    for alpha in [0,5,10,15,20,25]:
        res_q = run_xfoil('alfa', K_upper, K_lower, alpha, Re=Re, N_crit=N_crit, xtr_u=xtr_u, xtr_l=xtr_l, TE_gap=TE_gap)
        if res_q is not None:
            res_fastrun.append(res_q)
    
    cl_quick = [rs['cl'] for rs in res_fastrun]
    alpha_quick = [rs['alpha'] for rs in res_fastrun]

    clean_results = []
    if cl_design < max(cl_quick):
        slice_idx = [ i for i,vl in enumerate(cl_quick) if vl>cl_design ][0]
        if cl_quick[slice_idx] == cl_quick[slice_idx-1]:
            dist = 0.5
            alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
        else:
            dist = ( cl_design - cl_quick[slice_idx-1] )/( cl_quick[slice_idx] - cl_quick[slice_idx-1] )
            alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
        res_aguess = run_xfoil('alfa', K_upper, K_lower, alpha_design_guess, Re=Re, N_crit=N_crit, xtr_u=xtr_u, xtr_l=xtr_l, TE_gap=TE_gap)
        if abs(res_aguess['cl'] - cl_design)/cl_design <= 0.05:
            # pretty agressive tolerance of 5%, we should be in the linear region of the CL vs alpha curve
            res_075aguess = run_xfoil('alfa', K_upper, K_lower, 0.75*alpha_design_guess, Re=Re, N_crit=N_crit, xtr_u=xtr_u, xtr_l=xtr_l, TE_gap=TE_gap)
            clean_results.append(res_075aguess)
            clean_results.append(res_aguess)
            if all([cl_quick[iii+1] > cl_quick[iii] for iii in range(0,len(cl_quick)-1)]):
                # lift curve is still increasing, start at 10 deg
                alpha_start_search = 20.0
                res_stall_search_prev = res_a20
            else:
                # lift curve has rounded over, go to max-1 and search from there
                # max-1 is the last index we can guarentee is on the left of stall
                if cl_quick.index(max(cl_quick)) == 0:
                    alpha_start_search = 0.0
                    res_stall_search_prev = res_a00
                else:
                    alpha_start_search = alpha_quick[cl_quick.index(max(cl_quick))-1]
                    res_stall_search_prev = res_fastrun[cl_quick.index(max(cl_quick))-1]
            for iii in range(0,100):
                res_stall_search = run_xfoil('alfa', K_upper, K_lower, alpha_start_search+(iii+1)*0.5, Re=Re, N_crit=N_crit, xtr_u=xtr_u, xtr_l=xtr_l, TE_gap=TE_gap)
                clean_results.append(res_stall_search)
                if res_stall_search['cl'] > res_stall_search_prev['cl']:
                    # keep going
                    res_stall_search_prev = res_stall_search
                else:
                    # curve has rounded over
                    fastrunWorks = True
                    break

                if iii==99:
                    # could not find stall, revert to full run, though this will likely also fail
                    pass
        else:
            # no good, run full
            pass
    else:
        # could not find a good place to insert quick run marker, do full run
        pass

    if all([cr is not None for cr in clean_results]):
        # good
        pass
    else:
        fastrunWorks = False

    if fastrunWorks:
        return clean_results
    else:
        raise ValueError('Fastrun unsuccessful')

def airfoil_fitness(x):
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

        # design_matrix = [
        #     # tau,  CL,  spn,     Re
        #     [0.15, 1.2, 1.00, 10.0e6, ],
        #     [0.18, 1.2, 1.00, 10.0e6, ],
        #     [0.21, 1.2, 1.00, 12.0e6, ],
        #     [0.24, 1.2, 0.85, 13.0e6, ],
        #     [0.27, 1.2, 0.55, 16.0e6, ],
        # ]

        # design_matrix = [
        #     # tau,  CL,  spn,     Re
        #     [0.21, 1.4, 1.00, 12.0e6, ],
        # ]
        
        # ----------------------
        # unpack
        # ----------------------
        pid = x['pid']
        K_upper = x['individual'][0:int(len(x['individual'])/2)]
        K_lower = x['individual'][int(len(x['individual'])/2):]
        
        tau = x['tau']

        # ----------------------
        # reject self intersecting
        # ----------------------
        for i in range(0,len(K_upper)):
            if K_upper[i] <= K_lower[i]:
                return [pid, np.inf, False, -10] + [0]*15 + [0]*11
        
        # ----------------------
        # reject too large coefficients
        # ----------------------
        if max(abs(np.array(K_upper)))>2.0:
            return [pid, np.inf, False, -20] + [0]*15 + [0]*11
        if max(abs(np.array(K_lower)))>2.0:
            return [pid, np.inf, False, -30] + [0]*15 + [0]*11

        # ----------------------
        # build airfoil
        # ----------------------
        te_gap_lookup = {
            '15':  0.00196,
            '18':  0.00230,
            '21':  0.00262,
            '24':  0.00751,
            '27':  0.01012,
            '30':  0.01828,
            '33':  0.02644,
            '36':  0.02896,
        }

        afl_geo = Kulfan(TE_gap = te_gap_lookup[str(int(100*tau))])
        # afl_geo = Kulfan(TE_gap=0.0)
        afl_geo.upperCoefficients = K_upper
        afl_geo.lowerCoefficients = K_lower
        afl_geo.chord = 1.0*units.m
        # print(afl_geo.tau)
        # assert(abs(afl_geo.tau-tau)<0.01) # should never flag, but bad things happen if this is violated
        
        # ----------------------
        # Run Clean Data
        # ----------------------
        cl_design = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][1]
        Re = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][3]
        
        fastrun_clean_success = False
        try:
            clean_results = tryFastrun(K_upper, K_lower, cl_design, Re, N_crit=9.0, xtr_u=1.0, xtr_l=1.0, TE_gap = te_gap_lookup[str(int(100*tau))])
            fastrun_clean_success = True
        except:
            Re = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][3]
            lft = np.linspace(-30,-20,11)[0:-1]
            cnr = np.linspace(-20,20,81)
            rgt = np.linspace(20,30,11)[1:]
            alphaList = np.append(lft,np.append(cnr,rgt))
            
            clean_results=[]
            for alpha in alphaList:
                res = run_xfoil('alfa', K_upper, K_lower, alpha, Re=Re, N_crit=9.0, xtr_u=1.0, xtr_l=1.0, TE_gap = te_gap_lookup[str(int(100*tau))])
                if res is not None:
                    clean_results.append(res)

        cl_clean    = [en['cl'] for en in clean_results]
        cd_clean    = [en['cd'] for en in clean_results] 
        alpha_clean = [en['alpha'] for en in clean_results] 
        LoD_clean   = [en['cl']/en['cd'] for en in clean_results] 

        # ----------------------
        # Run Rough Data
        # ----------------------  
        try:
            if fastrun_clean_success:
                rough_results = tryFastrun(K_upper, K_lower, cl_design, Re, N_crit=3.0, xtr_u=0.05, xtr_l=0.05, TE_gap = te_gap_lookup[str(int(100*tau))])
            else:
                raise ValueError('Trivially trigger the except, clean quickrun was not successful')
        except:
            Re = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][3]
            lft = np.linspace(-30,-20,11)[0:-1]
            cnr = np.linspace(-20,20,81)
            rgt = np.linspace(20,30,11)[1:]
            alphaList = np.append(lft,np.append(cnr,rgt))
            
            rough_results=[]
            for alpha in alphaList:
                res = run_xfoil('alfa', K_upper, K_lower, alpha, Re=Re, N_crit=3.0, xtr_u=0.05, xtr_l=0.05, TE_gap = te_gap_lookup[str(int(100*tau))])
                if res is not None:
                    rough_results.append(res)

        cl_rough    = [en['cl'] for en in rough_results]
        cd_rough    = [en['cd'] for en in rough_results] 
        alpha_rough = [en['alpha'] for en in rough_results] 
        LoD_rough   = [en['cl']/en['cd'] for en in rough_results] 

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

            negative_peak_index_clean = mid_idx_clean
            for i in range(0, len(alpha_clean)):
                if cl_clean[mid_idx_clean-i] > cl_clean[mid_idx_clean-i-1]:
                    pass
                else:
                    negative_peak_index_clean = mid_idx_clean-i
                    break
                # will reach an index error if no peak is present

            negative_peak_index_rough = mid_idx_rough
            for i in range(0, len(alpha_rough)):
                if cl_rough[mid_idx_rough-i] > cl_rough[mid_idx_rough-i-1]:
                    pass
                else:
                    negative_peak_index_rough = mid_idx_rough-i
                    break
                # will reach an index error if no peak is present
        except:
            # threw an index error, so no peak is present
            return [pid, np.inf, False, -70] + [0]*15 + [0]*11
            
        # ----------------------
        # Find alpha_design
        # ----------------------
        cl_design = design_matrix[[dmr[0] for dmr in design_matrix].index(tau)][1]

        if cl_clean[positive_peak_index_clean] > cl_design:
            alpha_design = np.interp(cl_design, 
                                     np.array(cl_clean)[negative_peak_index_clean:positive_peak_index_clean], 
                                     np.array(alpha_clean)[negative_peak_index_clean:positive_peak_index_clean] )
        else:
            # raise ValueError("airfoil cannot reach target CL")
            return [pid, np.inf, False, -80] + [0]*15 + [0]*11
        
        # ----------------------
        # Begin objective function and constraints
        # ----------------------
        obj  = 0.0
        cons = []

        # ----------------------
        # Pure L/D performance
        # ----------------------
        LoD_clean_at_design_cl = np.interp(alpha_design, 
                                     np.array(alpha_clean)[negative_peak_index_clean:positive_peak_index_clean], 
                                     np.array(LoD_clean)[negative_peak_index_clean:positive_peak_index_clean] )

        LoD_rough_at_design_cl = np.interp(alpha_design, 
                                     np.array(alpha_rough)[negative_peak_index_rough:positive_peak_index_rough], 
                                     np.array(LoD_rough)[negative_peak_index_rough:positive_peak_index_rough] )
        
        # needed to shift negative L/D to be worse than positive
        obj +=  2.0*800/(500+LoD_clean_at_design_cl)
        obj += 24.0*800/(500+LoD_rough_at_design_cl)
        # equal weigting, baseline objective
        
        obj += 0.2*800/(500+max(LoD_clean))
        obj += 2.4*800/(500+max(LoD_rough))

        # ----------------------
        # Stall margin
        # ----------------------
        stall_margin_clean = alpha_clean[positive_peak_index_clean] - alpha_design
        stall_margin_rough = alpha_rough[positive_peak_index_rough] - alpha_design
        
        sm = 4.0
        # 0 if valid, otherwise 3.0-sm
        cons.append(sm-(min([sm,stall_margin_clean])))
        cons.append(sm-(min([sm,stall_margin_rough])))
        
        obj += (sm-min([sm, stall_margin_clean])) * 1e2
        obj += (sm-min([sm, stall_margin_rough])) * 1e2
        
        # ----------------------
        # lift margin
        # ----------------------
        lift_margin_clean = cl_clean[positive_peak_index_clean] - cl_design
        
        if lift_margin_clean >=0:
            # expected behavior
            obj += 0.5*lift_margin_clean
        else:
            # in some weird near-stall location, penalize
            obj += 0.5*abs(lift_margin_clean)*1e2

        # ----------------------
        # change in cl at alpha design due to roughness
        # ----------------------
        delta_cl_clean_to_rough_at_alpha_design = cl_design - np.interp(alpha_design, 
                                                                        np.array(alpha_rough)[negative_peak_index_rough:positive_peak_index_rough], 
                                                                        np.array(cl_rough)[negative_peak_index_rough:positive_peak_index_rough] )

        percent_delta_cl_clean_to_rough_at_alpha_design = delta_cl_clean_to_rough_at_alpha_design/cl_design
        
        # penalize if greater than 10%
        # assert(percent_delta_cl_clean_to_rough_at_alpha_design >=0)
        obj += (max([0.10, abs(percent_delta_cl_clean_to_rough_at_alpha_design)])-0.10) * 1e3
        
        # # ----------------------
        # # change in L/D at alpha design due to roughness
        # # ----------------------
        # delta_LoD_clean_to_rough_at_alpha_design = LoD_clean_at_design_cl - LoD_rough_at_design_cl
        # percent_delta_LoD_clean_to_rough_at_alpha_design = delta_LoD_clean_to_rough_at_alpha_design/LoD_clean_at_design_cl
     
        # # penalize if greater than 15%
        # if percent_delta_LoD_clean_to_rough_at_alpha_design >=0:
        #     # this is what it is supposed to be...
        #     obj += (max([0.15, percent_delta_LoD_clean_to_rough_at_alpha_design])-0.15) * 1e2
        # else:
        #     # something weird happened, increase penalty
        #     obj += abs(percent_delta_LoD_clean_to_rough_at_alpha_design) * 1e4

        # ----------------------
        # clean L/D curve fall off
        # ----------------------
        LoD_clean_15percent_left = np.interp(0.85*alpha_design, 
                                             np.array(alpha_clean)[negative_peak_index_clean:positive_peak_index_clean], 
                                             np.array(LoD_clean)[negative_peak_index_clean:positive_peak_index_clean] )

        LoD_clean_15percent_right = np.interp(1.15*alpha_design, 
                                             np.array(alpha_clean)[negative_peak_index_clean:positive_peak_index_clean], 
                                             np.array(LoD_clean)[negative_peak_index_clean:positive_peak_index_clean] )

        percent_change_LoD_clean_left  = (LoD_clean_15percent_left-LoD_clean_at_design_cl)/LoD_clean_at_design_cl
        percent_change_LoD_clean_right = (LoD_clean_15percent_right-LoD_clean_at_design_cl)/LoD_clean_at_design_cl
     
        # penalize if greater than 15%
        # assert(percent_delta_LoD_clean_to_rough_at_alpha_design >=0)
        obj += (max([0.15, abs(percent_change_LoD_clean_left )])-0.15) * 50 #1e2
        obj += (max([0.15, abs(percent_change_LoD_clean_right)])-0.15) * 50 #1e2

        # ----------------------
        # rough L/D curve fall off
        # ----------------------
        LoD_rough_15percent_left = np.interp(0.85*alpha_design, 
                                             np.array(alpha_rough)[negative_peak_index_rough:positive_peak_index_rough], 
                                             np.array(LoD_rough)[negative_peak_index_rough:positive_peak_index_rough] )

        LoD_rough_15percent_right = np.interp(1.15*alpha_design, 
                                             np.array(alpha_rough)[negative_peak_index_rough:positive_peak_index_rough], 
                                             np.array(LoD_rough)[negative_peak_index_rough:positive_peak_index_rough] )

        percent_change_LoD_rough_left  = (LoD_rough_15percent_left-LoD_rough_at_design_cl)/LoD_rough_at_design_cl
        percent_change_LoD_rough_right = (LoD_rough_15percent_right-LoD_rough_at_design_cl)/LoD_rough_at_design_cl

        # penalize if greater than 15%
        # assert(percent_delta_LoD_clean_to_rough_at_alpha_design >=0)
        obj += (max([0.15, abs(percent_change_LoD_rough_left )])-0.15) * 50 #1e2
        obj += (max([0.15, abs(percent_change_LoD_rough_right)])-0.15) * 50 #1e2
       
        # ----------------------
        # structure surrogates
        # ----------------------
        Ixx = afl_geo.Ixx.magnitude
        Iyy = afl_geo.Iyy.magnitude
        Izz = afl_geo.Izz.magnitude
        A   = afl_geo.area.magnitude

        Ixx_con = [0.00011000, 0.00017438, 0.00027518, 0.00041096, 0.00058321, 0.00079640, 0.00105795, 0.00137822 ]
        Iyy_con = [0.00397999, 0.00436351, 0.00493714, 0.00561409, 0.00633417, 0.00706380, 0.00779600, 0.00855043 ]
        Izz_con = [0.00408809, 0.00454606, 0.00521632, 0.00602287, 0.00691323, 0.00785849, 0.00885328, 0.00991577 ]
        A_con   = [0.08700496, 0.09995900, 0.11477620, 0.13051205, 0.14660942, 0.16289864, 0.17959744, 0.19731100 ]
        
        Ixx_target = Ixx_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        Iyy_target = Iyy_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        Izz_target = Izz_con[[dmr[0] for dmr in design_matrix].index(tau)] 
        A_target   = A_con[[dmr[0] for dmr in design_matrix].index(tau)] 

        # 0 if valid, otherwise target-val
        cons.append(Ixx_target-(min([Ixx_target,Ixx])))
        obj += (Ixx_target-min([Ixx_target, Ixx])) * 1e5
        
        cons.append(Iyy_target-(min([Iyy_target,Iyy])))
        obj += (Iyy_target-min([Iyy_target, Iyy])) * 1e4
        
        cons.append(Izz_target-(min([Izz_target,Izz])))
        obj += (Izz_target-min([Izz_target, Izz])) * 1e4
        
        cons.append(A_target-(min([A_target,A])))
        obj += (A_target-min([A_target, A])) * 1e4
        
        # ----------------------
        # reject leading edges with too tight radii
        # ----------------------
        leading_edge_radius_upper, leading_edge_radius_lower = afl_geo.leadingEdgeRadius()
        ler_con = [0.007, 0.008, 0.01, 0.025, 0.03, 0.04, 0.06, 0.08]
        
        leruViolation = -1*min([0,leading_edge_radius_upper-ler_con[[dmr[0] for dmr in design_matrix].index(tau)]])
        lerlViolation = -1*min([0,leading_edge_radius_lower-ler_con[[dmr[0] for dmr in design_matrix].index(tau)]])
        
        cons.append(leruViolation)
        cons.append(lerlViolation)
        
        obj += leruViolation * 1e4
        obj += lerlViolation * 1e4

        # ----------------------
        # reject violations of TE_cone
        # ----------------------
        cone_angle = 10
        te_frac = 0.98
        
        height_upper_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaUpper)
        height_lower_at_98percent = np.interp(te_frac,afl_geo.psi,afl_geo.zetaLower)

        midpoint_at_98_percent = ( height_upper_at_98percent + height_lower_at_98percent )/2

        upper_10deg_cone = midpoint_at_98_percent + np.tan(cone_angle/2/180*np.pi)*(1-te_frac)
        lower_10deg_cone = midpoint_at_98_percent - np.tan(cone_angle/2/180*np.pi)*(1-te_frac)

        te_valid = True
        teViolation = 0
        for i,psi_val in enumerate(afl_geo.psi):
            if psi_val >= te_frac:
                zeta_upper_val = afl_geo.zetaUpper[i]
                zeta_lower_val = afl_geo.zetaLower[i]
                h_upper = upper_10deg_cone - (psi_val-te_frac)*(upper_10deg_cone/(1-te_frac))
                h_lower = lower_10deg_cone - (psi_val-te_frac)*(lower_10deg_cone/(1-te_frac))

                if zeta_upper_val < h_upper:
                    teViolation += (h_upper-zeta_upper_val)                
                if zeta_lower_val > h_lower:
                    teViolation += (zeta_lower_val-h_lower)
        
        cons.append(teViolation)
        obj += teViolation*1e4            
       

        # ----------------------
        # reject if the max thickness is further forward than 25%
        # ----------------------
        tau_loc = afl_geo.taumax_psi
        if tau_loc < 0.30:
            cons.append(0.30-tau_loc)
            obj += (0.30-tau_loc) * 1e3
        else:
            cons.append(0)

        # ----------------------
        # attempt to prevent a bulge in the lower aft section
        # ----------------------
        bsl_cutoff = 0.8
        zeta_05 =np.interp(bsl_cutoff, afl_geo.psi, afl_geo.zetaLower) 
        er_tot = 0
        for i,psv in enumerate(afl_geo.psi):
            if psv >= bsl_cutoff:
                er = afl_geo.zetaLower[i] - np.interp(psv, [bsl_cutoff, 1.0], [zeta_05, afl_geo.zetaLower[-1]])
                er_tot += abs(min([0,er]))
        obj += er_tot*1e3

        # ----------------------
        # attempt to stop significantly asymmetric leading edges
        # ----------------------
        radii = afl_geo.leadingEdgeRadius()
        er = abs(radii[0]-radii[1])/min(radii)
        if er >=0.5:
            obj += 2*er


        # ----------------------
        # fix te wedge
        # ----------------------
        te_gap_diff = afl_geo.upperCoefficients[-1] - afl_geo.lowerCoefficients[-1]
        if te_gap_diff < 0.05:
            obj += (0.05-te_gap_diff) *1e4
        if te_gap_diff > 0.18:
            obj += te_gap_diff *1e4


        # ----------------------
        # Verify thickness, though this should be true by construction (except when nans were introduced)
        # ----------------------
        # negate because 0 is an unviolated constraint
        cons.append(int(not abs(afl_geo.tau-tau)<1e-4))

        # ----------------------
        # return
        # ----------------------
        con_tag = all([c==0 for c in cons])
        
        return [
                pid,
                obj,
                con_tag,
                alpha_design,
                LoD_clean_at_design_cl,
                stall_margin_clean,
                stall_margin_rough,
                lift_margin_clean,
                delta_cl_clean_to_rough_at_alpha_design,
                None,# delta_LoD_clean_to_rough_at_alpha_design, # this has been removed
                LoD_clean_15percent_left,
                LoD_clean_15percent_right,
                afl_geo.tau.magnitude,
                leading_edge_radius_upper, 
                leading_edge_radius_lower,
                Ixx,
                Iyy,
                Izz,
                A,
            ] + cons
    except:
        return [pid, np.inf, False, -90] + [0]*15 + [0]*11
