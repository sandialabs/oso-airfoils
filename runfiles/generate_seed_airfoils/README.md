Generating the Seed Airfoils
----------------------------

To improve the quality of guesses and rate of convergence, we pre-sample for airfoil geometries that are within a reasonably close distance of feasible designs.  This folder as it appears here is fully post-processed and needs no modification.  However, to replicate the work the following steps may be followed.  This assumes the installation of `mpi4py` using the command `pip install mpi4py` or similar as appropriate.

1.  Run the `find_seeds.py` script using the command `mpirun -n 188 python -m mpi4py find_seeds.py` or similar.  This produces the `viable_airfoils` folder, which consists of a number of `.txt` files that contain the Kulfan coefficients of the found airfoils we use to seed our optimization runs
    - Note: this may be run on a smaller number of processors (eg. `mpirun -n 188 python -m mpi4py find_seeds.py`) or on a single thread (eg. `python find_seeds.py`), however this may take substantial time.
2. Running `collect_viable_airfoils.py` collects the viable airfoils into a `all_viable_airfoils.txt` file that is used in the `newMember` function of the optimization run.  This also produces two plots that may be of interest.