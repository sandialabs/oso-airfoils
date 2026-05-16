Running the Optimization
------------------------

This folder contains all of the relevant files necessary to create your own airfoils or to replicate the work performed in this project.

A typical workflow is as follows:

1. Choose or create an input file (JSON or YAML format; see the Input File section below):
   - `quickstart.json` — lightweight neuralfoil run (N_k=10, N_pop=200, tau=0.21); good starting point for a personal machine
   - `run_001.json`–`run_008.json` — full xfoil-based runs, one per thickness (tau=0.15–0.36)
   - `run_nf_001.json`–`run_nf_008.json` — full neuralfoil-based runs, one per thickness
   - `example.yaml` — heavily annotated reference file documenting every available parameter
2. Run the optimization by passing the input file to `common_runner.py`:
   - **xfoil** (HPC, ~188 cores): `mpirun -n 188 python -m mpi4py common_runner.py run_001.json`
   - **neuralfoil** (HPC, ~96 cores): `mpirun -n 96 python -m mpi4py common_runner.py run_nf_001.json`
   - **personal machine**: `mpirun -n 8 python -m mpi4py common_runner.py quickstart.json`
3. Monitor per-generation Pareto front output printed to stdout.  Population snapshots are saved every generation into a timestamped output folder.
4. To regenerate all pre-built JSON files and batch shell scripts, run `generate_jsons.py` (xfoil) or `generate_jsons_nf.py` (neuralfoil).  These scripts also regenerate `all_run.sh` and `runcases.txt`.
5. To batch-submit all thickness cases at once, use `all_run.sh` (xfoil) or `all_run_nf.sh` (neuralfoil).

A personal laptop or computer running neuralfoil with N_k=8 and N_pop=200 will finish in roughly 24–48 hours and produce a reasonable result.  Using xfoil requires considerably more cores (128+) to run in a practical time.


Output
------

Each run creates an output folder under `outfile_leader` (default: `./`) named:

```
c{case}_t{tau*100}_k{N_k}_n{N_pop}_l{CL*10}_e{Re/1e5}__{timestamp}/
```

For example: `c112_t21_k16_n752_l15_e120__2026_05_15_23-56-1809/`

Inside this folder, the runner saves:

- A copy of the input JSON/YAML file at the start of the run
- Per-generation population snapshots: `population_{filecode}_g{generation}.json`

Each population snapshot is a JSON file containing:
- `input_parameters` — the full parameter dictionary used for this run, including start time, the generation number at time of write, and the path to the runner script
- `population` — a list of population members, each containing:
  - `K_upper` — list of N_k/2 Kulfan upper-surface coefficients
  - `K_lower` — list of N_k/2 Kulfan lower-surface coefficients
  - Computed quantities per member (objectives, constraint violations, aerodynamic metrics): `obj1`, `obj2`, `con_tag`, `alpha_design`, `LoD_clean_at_design`, `LoD_rough_at_design`, `stall_margin_clean`, `stall_margin_rough`, `lift_margin_clean`, `delta_cl_from_roughness`, `Ixx`, `Iyy`, `Izz`, `A`, `cpmin`, constraint tags, and `pareto_index`

The `pareto_index` field equals 1 for Pareto-front members and higher integers for lower-ranked fronts.  Each generation, the runner prints a table of the top Pareto-front members sorted by clean L/D.


Input File Format
-----------------

`common_runner.py` accepts either a `.json` or `.yaml`/`.yml` input file.  Every parameter available is documented with inline comments in `example.yaml`; that file is the authoritative reference.  The sections below summarize the parameters by category.

**Required parameters:** `case_number`, `tau`, `N_k`, `N_pop`, `N_generations`, `CL`, `Re`, `tool`, and all constraint weightings.  The pre-generated JSON files already contain a complete, valid set of every required parameter.


Core Optimization Parameters
-----------------------------

| Parameter | Description |
| :-------- | :---------- |
| `case_number` | Integer identifier used in output filenames |
| `tau` | Thickness-to-chord ratio; supported values: 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36 |
| `N_k` | Total number of Kulfan design variables, split evenly between upper and lower surfaces.  N_k=8 gives 4 coefficients per surface.  Minimum of 4 (2+2); recommend at least 8; limited gains above 16.  Higher values allow more shape flexibility but can over-tailor to the design point. |
| `N_pop` | GA population size.  Recommend at least 20 × N_k (e.g., ~400 for N_k=16, ~200 for N_k=8).  Must be divisible by 4. |
| `N_generations` | Number of GA generations to run before terminating.  Recommend ~1000 for production runs. |
| `CL` | Design-point lift coefficient (see defaults table below) |
| `Re` | Design-point Reynolds number (see defaults table below) |
| `tool` | Aerodynamic solver: `"xfoil"` or `"neuralfoil"` |
| `outfile_leader` | Directory where the output folder is created (default: `"./"`) |
| `xfoil_path` | Path to the xfoil executable (default: searches `$PATH`) |
| `xfoil_tempfile_path_leader` | Path prefix for xfoil temp files (default: `"t_"`) |


Aerodynamic Analysis Parameters
--------------------------------

Two polar sweeps are run for every airfoil candidate: a **clean** case representing attached laminar flow and a **rough** case representing leading-edge contamination.  Both use the same design Reynolds number.

| Parameter | Default | Description |
| :-------- | :-----: | :---------- |
| `N_crit_clean` | 9.0 | Transition amplification factor for clean analysis (free stream turbulence) |
| `xtp_u_clean` | 1.0 | Upper surface forced-transition location for clean case (1.0 = free transition) |
| `xtp_l_clean` | 1.0 | Lower surface forced-transition location for clean case |
| `alpha_min_clean` | 0 | Minimum angle of attack for clean sweep (degrees) |
| `alpha_max_clean` | 30 | Maximum angle of attack for clean sweep; must exceed stall |
| `alpha_step_clean` | 1 | Angle of attack step size for clean sweep (degrees) |
| `N_crit_rough` | 3.0 | Transition amplification factor for rough analysis (higher turbulence) |
| `xtp_u_rough` | 0.05 | Upper surface forced-transition location for rough case (5% chord) |
| `xtp_l_rough` | 0.05 | Lower surface forced-transition location for rough case (5% chord) |
| `alpha_min_rough` | 0 | Minimum angle of attack for rough sweep (degrees) |
| `alpha_max_rough` | 20 | Maximum angle of attack for rough sweep |
| `alpha_step_rough` | 1 | Angle of attack step size for rough sweep (degrees) |
| `xfoil_timelimit` | 15 | Time limit in seconds for a single xfoil case before it is terminated (ignored by neuralfoil) |
| `neuralfoil_model` | `"xxlarge"` | NeuralFoil model size; one of `xxsmall`, `xsmall`, `small`, `medium`, `large`, `xlarge`, `xxlarge`, `xxxlarge`.  Larger models are more accurate but slower. |
| `N_tries` | 1 | Number of times to retry a failed analysis before giving up; primarily useful for xfoil |
| `N_points_moi` | 20 | Number of airfoil surface points used to compute moments of inertia; 20 gives <1% error |


Objectives and Constraint Thresholds
--------------------------------------

The optimizer simultaneously maximizes **clean L/D** and **rough L/D** at the design angle of attack (the angle at which the clean polar achieves CL_design).  This is a two-objective Pareto problem solved with NSGA-II.

Constraint violations are added as a penalty to both objective values, collapsing the multi-objective problem into a single constraint-weighted penalized space.  A separate boolean `con_tag` tracks whether all constraints are simultaneously satisfied.

| Parameter | Default | Description |
| :-------- | :-----: | :---------- |
| `target_stall_margin` | 4.0 | Minimum required stall margin in degrees (alpha_stall − alpha_design), checked in both clean and rough conditions |
| `percent_delta_cl_from_roughness_threshold` | 0.10 | Maximum tolerated fractional drop in CL at alpha_design due to roughness (10%) |
| `percent_LoD_falloff_threshold` | 0.15 | Maximum tolerated fractional drop in L/D at alpha_design ± `alpha_falloff_offset`, evaluated in both clean and rough conditions (4 terms total) |
| `alpha_falloff_offset` | 1.0 | Angle of attack offset in degrees used to assess L/D curve breadth |
| `cl_max_limit_clean` | null | Hard upper limit on peak clean CL (useful for stall-regulated turbines); null disables |
| `cl_max_limit_rough` | null | Hard upper limit on peak rough CL; null disables |
| `target_cl` | null | Secondary rough-condition CL floor at a specified alpha (used for thick airfoils); null disables |
| `target_alpha` | null | Alpha at which `target_cl` is enforced in the rough polar; must be set if `target_cl` is set |
| `CMc_min` | null | Maximum allowable clean pitching moment magnitude in the band alpha_design ± `cm_alpha_band`; null disables |
| `CMr_min` | null | Maximum allowable rough pitching moment magnitude; null disables |
| `cm_alpha_band` | 5.0 | Half-width of the alpha band (degrees) over which moment constraints are checked |
| `cp_min_design` | null | Minimum allowable Cp anywhere on the airfoil at alpha_design (clean and rough); useful for cavitation studies |
| `cp_min_at_alpha_offset` | null | Minimum allowable Cp at alpha_design + `cp_min_alpha_offset` |
| `cp_min_alpha_offset` | null | Alpha offset in degrees for the `cp_min_at_alpha_offset` constraint |
| `cp_min_prestall` | null | Minimum allowable Cp anywhere on the airfoil, evaluated at all pre-stall alphas |


Structural Constraints
-----------------------

Structural constraints are defined as minimum values of cross-section properties non-dimensionalized by chord.  The structural surrogate properties are Ixx (edgewise), Iyy (flatwise), Izz (polar), and enclosed area A.  If a constraint value is set to `null`, an internally fitted geometry function (`geometry_functions.py`) fills in a default based on tau, calibrated to match the OSO design family.

| Parameter | Description |
| :-------- | :---------- |
| `TE_gap` | Trailing edge gap as a fraction of chord (exact; used in geometry construction). Default from geometry function if null. |
| `cone_angle` | Half-included angle (degrees) of the trailing-edge keep-out cone, centered on the camber line and applied aft of `te_frac`.  Prevents excessively blunt or diverging trailing edges.  Recommended: 10° for tau ≤ 0.18, 5° for 0.18–0.27, 0° for tau ≥ 0.30. |
| `te_frac` | Chord fraction where the TE cone constraint begins to be enforced (default: 0.95) |
| `Ixx_con` | Minimum edgewise moment of inertia |
| `Iyy_con` | Minimum flatwise moment of inertia |
| `Izz_con` | Minimum polar moment of inertia |
| `A_con` | Minimum enclosed cross-section area |
| `ler_con` | Minimum leading edge radius (applied to both upper and lower) |
| `ler_skew_factor` | Maximum ratio of the larger to smaller leading edge radius (default: 1.9).  Prevents highly asymmetric leading edges. |
| `max_thickness_loc` | Most forward chord fraction where the overall maximum thickness may be located (default: 0.275) |
| `max_thickness_loc_upper` | Most forward chord fraction for the upper-surface maximum thickness (default: 0.275) |
| `max_thickness_loc_lower` | Most forward chord fraction for the lower-surface maximum thickness (default: 0.275) |
| `min_radius_location_upper` | Minimum chord fraction for the location of minimum radius of curvature on the upper surface, within the first `min_radius_location_cutoff` of chord.  Null applies a geometry-function default. |
| `min_radius_location_lower` | Same as above for lower surface |
| `min_radius_location_cutoff` | Chord fraction defining the near-leading-edge search region for the radius-of-curvature location constraint (default: 0.08) |
| `curvature_bound` | Upper limit on second derivative of the upper surface aft of `ec_cutoff` (default: −750).  Prevents highly concave upper surfaces near the trailing edge. |
| `ec_cutoff` | Chord fraction marking the start of the aft-curvature enforcement region (default: 0.90) |
| `toothpick_height` | Minimum non-dimensional height at `toothpick_location` (default: 0.01).  Prevents needle-like airfoils. |
| `toothpick_location` | Chord fraction where the toothpick constraint is applied (default: 0.85) |

The OSO values for each tau are:

| Tau | TE_gap  | Ixx_con    | Iyy_con    | Izz_con    | A_con      | ler_con | cone_angle |
| :-: | :-----: | :--------: | :--------: | :--------: | :--------: | :-----: | :--------: |
| 15  | 0.00196 | 0.00011000 | 0.00397999 | 0.00408809 | 0.08700496 | 0.007   | 10°        |
| 18  | 0.00230 | 0.00017438 | 0.00436351 | 0.00454606 | 0.09995900 | 0.008   | 10°        |
| 21  | 0.00262 | 0.00027518 | 0.00493714 | 0.00521632 | 0.11477620 | 0.010   |  5°        |
| 24  | 0.00751 | 0.00041096 | 0.00561409 | 0.00602287 | 0.13051205 | 0.025   |  5°        |
| 27  | 0.01012 | 0.00058321 | 0.00633417 | 0.00691323 | 0.14660942 | 0.030   |  5°        |
| 30  | 0.01140 | 0.00079640 | 0.00706380 | 0.00785849 | 0.16289864 | 0.040   |  0°        |
| 33  | 0.01140 | 0.00105795 | 0.00779600 | 0.00885328 | 0.17959744 | 0.060   |  0°        |
| 36  | 0.01140 | 0.00137822 | 0.00855043 | 0.00991577 | 0.19731100 | 0.080   |  0°        |


Aerodynamic Design Defaults
----------------------------

These are the CL and Re defaults applied when those keys are not set in the JSON file.

| Tau |  Re   | CL  |
| :-: | :---: | :-: |
| 15  | 12e6  | 1.5 |
| 18  | 12e6  | 1.5 |
| 21  | 12e6  | 1.5 |
| 24  | 13e6  | 1.4 |
| 27  | 16e6  | 1.3 |
| 30  | 18e6  | 1.2 |
| 33  | 16e6  | 1.2 |
| 36  | 13e6  | 1.2 |


Constraint Weightings
----------------------

Constraint violations are penalized by multiplying the violation magnitude by the corresponding weight and adding the result to both objective values.  Larger weights make a constraint harder to violate; setting a weight to `0.0` or `null` disables that constraint entirely.

The weights used in the pre-generated JSON files reflect the OSO design study and are a reasonable starting point.  The full list with descriptions is in `example.yaml`.  Notable entries:

| Parameter | Default | Notes |
| :-------- | :-----: | :---- |
| `stall_margin_clean_weighting` | 1e2 | Clean stall margin |
| `stall_margin_rough_weighting` | 1e2 | Rough stall margin |
| `lift_margin_clean_weighting` | 0.5 | Soft incentive to maximize CL_max above CL_design (not a hard constraint); set to 0 if `cl_max_limit_clean` is used |
| `delta_cl_from_roughness_weighting` | 1e4 | Roughness CL sensitivity |
| `LoD_falloff_weighting` | 50 | L/D breadth (applied to all 4 terms) |
| `ixx_weighting` | 1e6 | Edgewise inertia (weight is large due to small numerical magnitude of Ixx) |
| `iyy_weighting` | 1e4 | Flatwise inertia |
| `izz_weighting` | 1e4 | Polar inertia |
| `a_weighting` | 1e4 | Cross-section area |
| `leading_edge_radius_upper_weighting` | 1e3 | Upper LE radius |
| `leading_edge_radius_lower_weighting` | 1e3 | Lower LE radius |
| `te_cone_violation_weighting` | 1e5 | TE cone (large weight needed due to small violation magnitudes) |
| `curvature_weighting` | 100 | Upper surface concavity |
| `max_thickness_lower_weighting` | 5e4 | Lower surface thickness location (higher than upper due to strong optimizer pressure) |
| `infeasibility_penalty` | 1e4 | Blanket penalty added to any infeasible design |


Genetic Algorithm Parameters
------------------------------

| Parameter | Default | Description |
| :-------- | :-----: | :---------- |
| `probability_of_mutation` | 0.3 | Per-gene bit-flip probability.  The current implementation represents each design variable as a binary-encoded floating-point number; a value of 0.3 means 30% of bit positions are eligible to flip.  This is intentionally high to encourage diversity. |
| `N_mutations` | 4 | Number of gene mutations applied to a child if it is selected for mutation (all-or-nothing selection) |
| `N_crossovers` | 3 | Number of crossover points during chromosome recombination; recommend 3 |
| `maximum_parent_fraction` | 0.7 | Maximum fraction of the parent population that can survive into the next generation; limits elitism |
| `front1_cap_fraction` | 0.5 | Maximum fraction of the population allowed to come from Pareto front 1.  When front 1 grows too large, NSGA-II degenerates to crowding-distance selection and diversity collapses.  This cap forces excess front-1 members to be displaced by members from lower fronts. |

The NSGA-II implementation is in `ga_new_generation_mpi_nsga2_v2.py`.  Two fixes versus the original NSGA-II are applied: (1) duplicate design-vector removal from the combined parent+child pool before non-dominated sorting, and (2) the front-1 population cap described above.


Continuation Runs
-----------------

Any run can be resumed from a saved population snapshot by setting `continuation_file` in the input file:

```json
"continuation_file": "./c112_t21_k16_n752__2026_05_10/population_c112_t21_k16_n752_g300.json",
"continuation_file_overwrite": false,
"N_generations": 500
```

- `continuation_file` may be either a path to a specific population JSON file, or a path to the output folder (in which case the alphabetically last JSON file in that folder is used, typically the most recent generation).
- When `continuation_file_overwrite` is `false` (default), all optimization parameters are taken from the saved file and `N_generations` is interpreted as **additional** generations beyond those already run.
- When `continuation_file_overwrite` is `true`, the optimization parameters in the current input file overwrite those from the saved file.  The population is still loaded from the saved file.
- Output continues writing to the same output folder.


Generating Input Files Programmatically
-----------------------------------------

`generate_jsons.py` and `generate_jsons_nf.py` produce a full set of per-thickness input files and the corresponding `all_run.sh` batch script.  Edit the `default_dict` near the top of either file to change shared parameters, then run the script.  The `tau_data_dict` within each script sets the per-thickness structural constraints, trailing edge gap, and aerodynamic design point.

The distinction between the two scripts is the `tool` field: `generate_jsons.py` sets `"tool": "xfoil"` and `generate_jsons_nf.py` sets `"tool": "neuralfoil"`.


Modifying the Objective Function
----------------------------------

The objective function and all constraint logic is in `wt_objective_nsga2.py`.  The MPI task distribution and GA loop are in `common_runner.py` and `ga_new_generation_mpi_nsga2_v2.py`.  The initial seed population is drawn from a library of 100 pre-generated airfoils in `newMember.py`, scaled to the target tau before the first generation.
