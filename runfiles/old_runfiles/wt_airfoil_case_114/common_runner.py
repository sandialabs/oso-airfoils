# mpirun -n 8 python -m mpi4py c63_t15_l15_r124_k8_n200.py

import math
import shutil
from wt_objective_nsga2 import airfoil_fitness
# from ga_new_generation_mpi_nsga2 import newGeneration
from ga_new_generation_mpi_nsga2_v2 import newGeneration
from kulfan import Kulfan
import numpy as np
import time
import os
import sys
import copy
import json
import yaml
import argparse
from newMember import newMember
import platform
import subprocess

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

def save_json(fname, pop, params, labels, datestr, current_generation):
    save_dict = {}
    save_dict['input_parameters'] = params
    save_dict['input_parameters']['start_time'] = datestr
    save_dict['input_parameters']['current_generation'] = current_generation
    save_dict['input_parameters']['write_time'] = time.strftime("%Y_%m_%d_%H-%M-%S", time.localtime()) + f"{int((time.time() % 1) * 100):02d}"
    save_dict['input_parameters']['path_to_here'] = str(path_to_here)
    save_dict['input_parameters']['operating_system'] = platform.system()
    pop_arr = np.array(pop)

    # Split the per-row label list into upper-coeff, lower-coeff, and the rest.
    n_half = sum(1 for lb in labels if lb.startswith('U') and lb[1:].isdigit())
    rest_labels = labels[2 * n_half:]

    population = []
    for row in pop_arr:
        row = row.tolist()
        entry = {
            'K_upper': row[:n_half],
            'K_lower': row[n_half:2 * n_half],
        }
        entry.update(dict(zip(rest_labels, row[2 * n_half:])))
        population.append(entry)

    save_dict['population'] = population
    with open(fname, 'w') as f:
        json.dump(save_dict, f, indent=4)

    return save_dict
    # np.savetxt(ldr + folderstr + '/population_%s_g%d.txt'%(filecode,i),np.array(pop))

# =========================================================================================================
# =========================================================================================================

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run airfoil optimization with JSON or YAML input')
parser.add_argument('input_file', help='JSON or YAML input file containing parameters')
args = parser.parse_args()

# Load parameters from JSON or YAML file based on extension
def read_input_file(input_file):
    _ext = os.path.splitext(input_file)[1].lower()
    with open(input_file, 'r') as f:
        if _ext in ('.yaml', '.yml'):
            return yaml.safe_load(f)
        elif _ext == '.json':
            return json.load(f)
        else:
            raise ValueError(
                f"Unsupported input file extension '{_ext}'. Use .json, .yaml, or .yml."
            )
params = read_input_file(args.input_file)

# Check for required parameters
# required_params = ['case_number', 'tau', 'N_k', 'N_pop']
# missing_params = [param for param in required_params if param not in params or params[param] is None]
# if missing_params:
#     raise ValueError(f"Missing required parameters in JSON file: {missing_params}")

if 'continuation_file_overwrite' not in params:
    params['continuation_file_overwrite'] = False
if params['continuation_file_overwrite'] is None:
    params['continuation_file_overwrite'] = False

is_continuation = False
if 'continuation_file' in params:
    if params['continuation_file'] is not None:
        if params.get('continuation_file_overwrite') is False:
            if os.path.isdir(params['continuation_file']):
                # list all json files in the directory and pick the last one when sorted
                json_files = [f for f in os.listdir(params['continuation_file']) if f.endswith('.json')]
                if not json_files:
                    raise ValueError(f"No JSON files found in directory specified by 'continuation_file': {params['continuation_file']}")
                json_files.sort()
                continuation_file = os.path.join(params['continuation_file'], json_files[-1])
            else:
                continuation_file = params['continuation_file']

            assert(os.path.isfile(continuation_file))

            is_continuation = True
            params_original = copy.deepcopy(params)
            confile_data = read_input_file(params_original['continuation_file'])
            params = confile_data['input_parameters']

# Extract parameters with defaults
case_number = params.get('case_number')
tau = params.get('tau')
N_k = params.get('N_k')
N_pop = params.get('N_pop')
CL = params.get('CL')
Re = params.get('Re')
file_system = params.get('file_system')
N_generations = int(params.get('N_generations'))
ldr = params.get('outfile_leader')
te_gap = params.get('TE_gap')

# Validate file_system parameter
if file_system is not None and file_system not in [0,1,2,3]:
    raise ValueError('filesystem flag must be 0 (default), 1 (gpfs), 2 (pscratch), or 3 (tscratch)')

# Create filecode from parameters for output naming
filecode = f"c{case_number}_t{int(tau*100)}_k{N_k}_n{N_pop}"
if CL is not None:
    filecode += f"_l{int(CL*10)}"
if Re is not None:
    filecode += f"_e{int(Re/1e5)}"
if file_system is not None:
    filecode += f"_s{file_system}"

if ldr is None:
    ldr = '.'+os.sep

if N_k < 4:
    raise ValueError("Must use at least 2 design variables top and bottom")

previous_generation_count = 0
if is_continuation:
    previous_generation_count = int(params.get('current_generation'))
    N_generations = previous_generation_count + int(params_original.get('N_generations'))
# =========================================================================================================
# =========================================================================================================

labels = []

for i in range(0,int(N_k/2)):
    labels += ['U%d'%(i+1)]

for i in range(0,int(N_k/2)):
    labels += ['L%d'%(i+1)]

labels += [
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
        'pareto_index',
    ]




    

# =========================================================================================================
# =========================================================================================================

pop = None 

if rank == 0:
    if is_continuation:
        folderstr = params_original['continuation_file'].split(os.sep)[-2]
        datestr = params['start_time']
    else:
        datestr = time.strftime("%Y_%m_%d_%H-%M-%S", time.localtime()) + f"{int((time.time() % 1) * 100):02d}"
        folderstr = filecode + '__' + datestr
        if not os.path.isdir(ldr + folderstr):
            os.mkdir(ldr + folderstr) 
        shutil.copy(args.input_file, ldr + folderstr + os.sep + datestr + '.' + args.input_file.split('.')[-1])
        shutil.move(args.input_file, ldr + folderstr + os.sep + args.input_file)

if rank == 0:
    if is_continuation:
        # Reconstruct the population array from the saved JSON. Each entry has
        # 'K_upper' and 'K_lower' lists followed by the remaining scalar labels;
        # save_json wrote rows as: K_upper + K_lower + [entry[lbl] for lbl in rest_labels].
        rest_labels = labels[N_k:]
        pop = []
        for entry in confile_data['population']:
            row = list(entry['K_upper']) + list(entry['K_lower'])
            row += [entry[lbl] for lbl in rest_labels]
            pop.append(row)
        pop = np.array(pop)
    else:
        pop = newMember(int(N_k/2),tau,N_pop,te_gap = te_gap)

pop = comm.bcast(pop, root=0)

if not is_continuation:
    pop = newGeneration(airfoil_fitness, pop, normalizationVector = [1]*N_k, encodingTypes=[float]*N_k, lowerBounds=[-2.0]*N_k, upperBounds=[2.0]*N_k, initalize=True, comm=comm, params=params)
    pop = comm.bcast(pop, root=0)

    if rank == 0:
        # np.savetxt(ldr + folderstr + '/population_%s_g0.txt'%(filecode),np.array(pop))
        save_json(ldr + folderstr + os.sep +'population_%s_g%s.json'%(filecode,str(0).zfill(math.ceil(np.log10(N_generations)))),np.array(pop),params,labels,datestr,0)

# =========================================================================================================
# =========================================================================================================
# if rank == 0:
counter = previous_generation_count+1

for i in range(previous_generation_count+1,int(1.2*N_generations)):
    pop_cache = None
    if rank == 0:
        cprint('Generation %d'%(counter))
        # Cache the previous population in the event of a crash
        pop_cache = copy.deepcopy(pop)
    pop_cache = comm.bcast(pop_cache, root=0)
    try:
        pop = newGeneration(airfoil_fitness, pop, normalizationVector = [1]*N_k, encodingTypes=[float]*N_k, lowerBounds=[-2.0]*N_k, upperBounds=[2.0]*N_k, initalize=False, comm=comm, params=params)
        pop = comm.bcast(pop, root=0)
    except:
        cprint('Error occurred during generation %d, reverting population to previous generation\n'%(counter))
        pop = pop_cache
        pop = comm.bcast(pop, root=0)
        continue
    if rank == 0:
        if len(pop) == 0:
           # Crash occurred, reload the previous generation and retry
           pop = pop_cache
           pop = comm.bcast(pop, root=0)
        # Save to file
        save_dict = save_json(ldr + folderstr + os.sep + 'population_%s_g%s.json'%(filecode,str(counter).zfill(math.ceil(np.log10(N_generations)))),np.array(pop),params,labels,datestr,counter)
        counter += 1

        pareto_elements = [elem for elem in save_dict['population'] if elem['pareto_index'] == 1]
        pareto_elements.sort(key=lambda elem: elem['LoD_clean_at_design'], reverse=True)
        tbl = {'Index':[], 'Clean L/D':[], 'Rough L/D':[], 'Feasible':[]}
        for ii, elem in enumerate(pareto_elements):
            tbl['Index'].append(ii)
            tbl['Clean L/D'].append(elem['LoD_clean_at_design'])
            tbl['Rough L/D'].append(elem['LoD_rough_at_design'])
            tbl['Feasible'].append(elem['con_tag'])
        pstr = ''
        for key in tbl.keys():
            pstr += key.ljust(15)
        pstr += '\n'
        n_pareto = len(pareto_elements)
        max_rows = 15
        if n_pareto <= max_rows:
            row_indices = list(range(n_pareto))
        else:
            row_indices = [int(round(j * (n_pareto - 1) / (max_rows - 1))) for j in range(max_rows)]
        for ii in row_indices:
            for key in tbl.keys():
                val = tbl[key][ii]
                if isinstance(val, float):
                    pstr += f"{val:.2f}".ljust(15)
                else:
                    pstr += str(val).ljust(15)
            pstr += '\n'

        cprint(pstr)

    should_break = comm.bcast(bool(rank == 0 and counter > N_generations), root=0)
    if should_break:
        break

        # np.savetxt(ldr + folderstr + '/population_%s_g%d.txt'%(filecode,i),np.array(pop))


