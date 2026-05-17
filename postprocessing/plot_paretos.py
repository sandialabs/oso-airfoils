import natsort
import numpy as np
import os

# Configuration for different constraint types
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

def get_default_folder_config():
    """Returns the default configuration of folder lists for different constraint types"""
    return {
        'non_constrained_folders': [
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_109/c109_t18_l15_k16_g2000_n752__2025_12_31_16-57',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_109/c109_t21_l15_k16_g2000_n752__2026_01_07_12-47',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_109/c109_t24_l14_k16_g2000_n752__2026_01_10_14-28',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_109/c109_t27_l13_k16_g2000_n752__2026_01_14_19-33',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_107/c107_t30_l12_k16_g2000_n752__2025_12_19_00-28',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_107/c107_t33_l12_k16_g2000_n752__2025_12_22_20-35',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_107/c107_t36_l12_k16_g2000_n752__2025_12_20_16-46',
        ],
        'moment_constrained_folders': [
            # '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_109/c109_t18_l15_k16_g2000_n752__2025_12_31_16-57',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t21_l15_k16_g2000_n752_m14_p14__2026_01_17_08-45',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t24_l14_k16_g2000_n752_m14_p14__2026_01_22_11-18',
            # "/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/c111_t24_l14_k16_g2000_n752__2026_03_07_14-35",
            # '/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_x3_s1__2026_02_26_09-48',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t27_l13_k16_g2000_n752_m14_p14__2026_02_21_12-43',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t30_l12_k16_g2000_n752_m14_p14__2026_02_27_13-34',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t33_l12_k16_g2000_n752_m14_p14__2026_03_05_10-37',
            # '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_110/c110_t36_l12_k16_g2000_n752_m14_p14__2026_03_11_13-39',
            '/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/c111_t36_l12_k16_g2000_n752_m14_p14__2026_03_12_23-29',
        ],
        'supercomputer_folders': [
            # "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t15_l15_k16_g2000_n752_x3_s1__2026_02_26_09-48",
            # "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t18_l15_k16_g2000_n752_x3_s2__2026_02_26_09-46",
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t18_l15_k16_g2000_n752_x3_s2__2026_03_09_09-55",
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t21_l15_k16_g2000_n752_x3_s3__2026_02_26_09-47",
            # "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_x3_s1__2026_02_26_09-48",
            '/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_x3_s1__2026_03_03_16-50/',
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t27_l13_k16_g2000_n752_x3_s2__2026_02_26_09-49",
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t30_l12_k16_g2000_n752_x3_s3__2026_02_26_09-54",
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t33_l12_k16_g2000_n752_x3_s1__2026_02_26_10-01",
            "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_101_to_110/case_110/c110_t36_l12_k16_g2000_n752_x3_s2__2026_02_26_10-12",
        ],
        'neuralfoil_folders': [
            "/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/c111_t24_l14_k16_g2000_n752__2026_03_07_14-35",
        ],
        'neuralfoil_postcompute': [
            "/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/population_c111_t24_l14_k16_g2000_n752_g882__2026_03_16_15-47"
            # "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_111/c111_t21_l15_k16_g2000_n752_x3_s3__2026_04_06_21-52"
            # "/Users/codykarcher/Dropbox/research/oso-airfoils/postprocessing/cases/cases_111_to_120/case_111/c111_t21_l11_k16_g2000_n752_x3_s2__2026_04_06_21-52"
        ],
    }

def process_airfoil_data(data_folders, constraint_type='non_constrained'):
    """Process airfoil data for either constrained or non-constrained cases"""
    
    family_data = {
        'tau' :      np.array([]),
        'psi':       np.array([]),
        'factor':    np.array([]),
        'cl':        np.array([]),
        're':        np.array([]),
        'camber':    np.array([]),
        'thickness': np.array([]),
        'rough':     np.array([]),
        'clean':     np.array([]),
    }

    for path_to_data in data_folders:
        # Extract tau value from path
        if 't36'in path_to_data:
            tau = '36'
        elif 't33'in path_to_data:
            tau = '33'
        elif 't30'in path_to_data:  
            tau = '30'
        elif 't27'in path_to_data:  
            tau = '27'
        elif 't24'in path_to_data:  
            tau = '24'
        elif 't21'in path_to_data:  
            tau = '21'
        elif 't18'in path_to_data:  
            tau = '18'
        elif 't15'in path_to_data:  
            tau = '15'
        else:
            continue

        # Set objective indices based on constraint type
        if constraint_type == 'non_constrained':
            if tau in ['18','21','24','27']:
                obj1ix = 22
                obj2ix = 21
            else:
                obj1ix = 21
                obj2ix = 20
        elif constraint_type == 'moment_constrained':
            if tau in ['18']:
                obj1ix = 22
                obj2ix = 21
            else:
                obj1ix = 21
                obj2ix = 20
        elif constraint_type == 'supercomputer':
            # if tau in ['18']:
            #     obj1ix = 22
            #     obj2ix = 21
            # else:
            obj1ix = 21
            obj2ix = 20
        elif 'neuralfoil' in constraint_type:
            obj1ix = 21
            obj2ix = 20
        else:
            raise ValueError("Invalid constraint type specified.")

        design_matrix = {
            # tau,  CL,  spn,     Re
            '15':  [ 1.5, 1.00, 10.0e6, ],
            '18':  [ 1.5, 1.00, 10.0e6, ],
            '21':  [ 1.5, 1.00, 12.0e6, ],
            '24':  [ 1.4, 0.85, 13.0e6, ],
            '27':  [ 1.3, 0.55, 16.0e6, ],
            '30':  [ 1.2, 0.50, 18.0e6, ],
            '33':  [ 1.2, 0.35, 16.0e6, ],
            '36':  [ 1.2, 0.20, 13.0e6, ],
        }

        if '_l' in path_to_data:
            l_ix_start = path_to_data.index('_l') + 2
            l_ix_end = path_to_data.index('_k')
            l_val = path_to_data[l_ix_start:l_ix_end]
            CL = float(l_val) * 0.1
        else:
            CL = design_matrix[tau][1]

        Re = design_matrix[tau][2]

        files = natsort.natsorted([f for f in os.listdir(path_to_data) if '.txt' in f], alg=natsort.ns.IGNORECASE)
        f = files[-1]
        data = np.loadtxt( os.path.join(path_to_data, f) )
        pareto_points = data[ data[:,-1] == 1 ]
        pareto_points = pareto_points[ np.argsort( pareto_points[:, obj1ix] ) ]
        rgh = pareto_points[:, obj1ix]

        for i in range(0, len(pareto_points)):
            if abs(pareto_points[i,16]) <= 0.1:
                afl = Kulfan( TE_gap=pareto_points[i,16] )
            else:
                afl = Kulfan( TE_gap=te_gap_lookup[tau] )
            afl.upperCoefficients = pareto_points[i, 0:8]
            afl.lowerCoefficients = pareto_points[i, 8:16]
            cb = (afl.zetaUpper + afl.zetaLower) / 2
            tk = afl.zetaUpper - afl.zetaLower
            psi = afl.psi
            fctr = (rgh[i] - min(rgh)) / (max(rgh) - min(rgh)) * np.ones( len(psi) )

            family_data['tau'] = np.concatenate( (family_data['tau'],  int(tau)/ 100 * np.ones( len(psi) ) ) )
            family_data['psi'] = np.concatenate( (family_data['psi'],  psi ) )
            family_data['factor'] = np.concatenate( (family_data['factor'],  fctr ) )
            family_data['cl'] = np.concatenate( (family_data['cl'],  CL * np.ones( len(psi) ) ) )
            family_data['re'] = np.concatenate( (family_data['re'],  Re * np.ones( len(psi) ) ) )
            family_data['camber'] =  np.concatenate( (family_data['camber'],  cb ) )
            family_data['thickness'] = np.concatenate( (family_data['thickness'],  tk ) )
            family_data['rough'] = np.concatenate( (family_data['rough'],  rgh[i] * np.ones( len(psi) ) ) )
            family_data['clean'] = np.concatenate( (family_data['clean'],  pareto_points[i, obj2ix] * np.ones( len(psi) ) ) )

    return family_data

def main(folder_config=None):
    """Main processing function that takes a dictionary of folder lists"""
    if folder_config is None:
        folder_config = get_default_folder_config()
    
    # Import required modules
    import pandas as pd
    import matplotlib.pyplot as plt
    from kulfan import Kulfan
    
    # Process datasets based on provided configuration
    processed_data = {}
    dataframes = {}
    
    # Define constraint type mapping (remove '_folders' suffix for constraint type)
    constraint_type_mapping = {
        'non_constrained_folders': 'non_constrained',
        'moment_constrained_folders': 'moment_constrained', 
        'supercomputer_folders': 'supercomputer',
        'neuralfoil_folders': 'neuralfoil',
        'neuralfoil_postcompute': 'neuralfoil_postcompute',
    }
    
    for folder_type, folders in folder_config.items():
        if folders:  # Only process if folder list is not empty
            constraint_type = constraint_type_mapping.get(folder_type, 'non_constrained')
            print(f"Processing {folder_type.replace('_folders', '')} data...")
            processed_data[folder_type] = process_airfoil_data(folders, constraint_type)
            dataframes[folder_type] = pd.DataFrame(processed_data[folder_type])
    
    # Save individual datasets
    filename_mapping = {
        'non_constrained_folders': 'airfoil_family_data.csv',
        'moment_constrained_folders': 'airfoil_family_data_mcon.csv',
        'supercomputer_folders': 'airfoil_family_data_supercomputer.csv',
        'neuralfoil_folders': 'airfoil_family_data_neuralfoil.csv',
        'neuralfoil_postcompute': 'airfoil_family_data_neuralfoil_postcompute.csv'
    }
    
    for folder_type, df in dataframes.items():
        filename = filename_mapping.get(folder_type, f'airfoil_family_data_{folder_type}.csv')
        df.to_csv(filename, index=False)
        print(f"Saved {filename} with {len(df)} points")
    
    return dataframes

def create_combined_plot(dataframes):
    """Create combined plot from processed dataframes"""
    import matplotlib.pyplot as plt
    
    # ===================================================================================================================
    # Combined Plot
    plt.figure(figsize=(14, 10))
    
    # Process each dataset type
    psi0_data = {}
    taus_by_type = {}
    
    for folder_type, df in dataframes.items():
        if len(df) > 0:
            psi0_data[folder_type] = df[df['psi'] == 0]
            taus_by_type[folder_type] = sorted(list(set(psi0_data[folder_type]['tau'].tolist())))

    # Get all unique taus for consistent coloring
    all_taus_lists = [taus for taus in taus_by_type.values() if taus]
    all_taus = sorted(list(set([tau for taus in all_taus_lists for tau in taus])))
    # colors = plt.cm.viridis(np.linspace(0, 1, len(all_taus)))
    colors = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa', '#ede13f', '#56b4ff', '#fca7c7', '#5d5d5d', '#000000']
    
    # Define plot styles for different dataset types
    plot_styles = {
        'non_constrained_folders': {'linestyle': '-', 'label_suffix': 'Non-constrained'},
        'moment_constrained_folders': {'linestyle': '--', 'label_suffix': 'Moment-constrained'},
        'supercomputer_folders': {'linestyle': ':', 'label_suffix': 'Supercomputer'},
        'neuralfoil_folders': {'linestyle': ':', 'label_suffix': 'Neuralfoil'},
        'neuralfoil_postcompute': {'linestyle': '-', 'label_suffix': 'Neuralfoil Postcompute'}
    }
    
    # Plot each dataset type
    for folder_type, tau_list in taus_by_type.items():
        if folder_type in psi0_data and len(psi0_data[folder_type]) > 0:
            style = plot_styles.get(folder_type, {'linestyle': '-', 'label_suffix': folder_type.replace('_folders', '')})
            
            for tau in all_taus:
                if tau in tau_list:
                    tau_data = psi0_data[folder_type][psi0_data[folder_type]['tau'] == tau]
                    color_idx = all_taus.index(tau)
                    if 'neuralfoil' not in folder_type:
                        plt.plot(tau_data['rough'], tau_data['clean'], style['linestyle'], 
                                color=colors[color_idx], linewidth=2.5, markersize=5,
                                label=f'τ = {tau:.2f} ({style["label_suffix"]})', alpha=0.8)
                    else:
                        plt.plot(tau_data['rough'], tau_data['clean'], style['linestyle'], 
                            color=[0,1,0], linewidth=2.5, markersize=5,
                            label=f'τ = {tau:.2f} ({style["label_suffix"]})', alpha=0.8)

    # Add styling
    plt.xlabel('Rough L/D', fontsize=14, fontweight='bold')
    plt.ylabel('Clean L/D', fontsize=14, fontweight='bold')
    plt.title('Airfoil Performance: Clean vs Rough L/D by Thickness Ratio\n(Solid: Non-constrained, Dashed: Moment-constrained, Dotted: Supercomputer)', 
              fontsize=16, fontweight='bold', pad=20)

    plt.legend(frameon=True, fancybox=True, shadow=True, 
              fontsize=10, loc='upper left', ncol=2)

    plt.grid(True, alpha=0.3, linestyle='--')

    # Add WT2 reference data
    WT2_data = {
        # This is RFOIL data, so not exactly right
        '21': (248.2, 120.2),  
        '24': (183.9, 118.1),  
        '27': (172.5, 112.2),  
        '30': (167.4, 97.4),   
        '33': (163.1, 87.1),   
        '36': (157.1, 75.4),   
    }

    WT2_data_xfoil = {
        '21': (255.8, 122.2),  
        '24': (175.2, 118.0),  
        '27': (182.0, 109.0),  
        '30': (177.8, 98.0),   
        '33': (166.6, 90.0),   
        '36': (159.7, 78.1),   
    }

    for tau, (clean_ld, rough_ld) in WT2_data.items():
        if float(tau)/100 in all_taus:
            color_idx = all_taus.index(float(tau)/100)
            # plt.plot(rough_ld, clean_ld, 'X', color=colors[color_idx], markersize=10, markeredgecolor='k', label=f'WT2 τ={tau}')
            plt.plot(rough_ld, clean_ld, 'o', color=colors[color_idx], markersize=10, markeredgecolor='k', label=f'WT2 τ={tau}')

    for tau, (clean_ld, rough_ld) in WT2_data_xfoil.items():
        if float(tau)/100 in all_taus:
            color_idx = all_taus.index(float(tau)/100)
            plt.plot(rough_ld, clean_ld, 'X', color=colors[color_idx], markersize=10, markeredgecolor='k', label=f'WT2 τ={tau}')
            # plt.plot(rough_ld, clean_ld, 'o', color=colors[color_idx], markersize=10, markeredgecolor='k', label=f'WT2 τ={tau}')

    plt.tight_layout()
    # plt.xlim([0, 150])
    # plt.ylim([0,400])

    
    # Style the axes
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_linewidth(1.5)
    plt.gca().spines['bottom'].set_linewidth(1.5)
    
    # plt.show()
    return plt.gcf()


# if __name__ == "__main__":
# Run with default configuration
from kulfan import Kulfan
import matplotlib.pyplot as plt

dataframes = main()
fig = create_combined_plot(dataframes)

print("Data processing and plotting complete!")
for folder_type, df in dataframes.items():
    print(f"{folder_type.replace('_folders', '')} data: {len(df)} points")