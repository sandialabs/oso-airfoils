Postprocessing
--------------

This folder contains all of the raw results of optimization runs, located in the `cases` folder.  Postprocessing is still being worked on.

The following scripts are included here, and are only valid for use with Cases 112 and greater:
- `compare_airfoils.py`--Contains a number of useful plotting features, most notably is the `compare_airfoils` function that allows for comparisons of airfoil performance
- `find_shape_functions.py`--Used to determine modified Chebyshev polynomial coefficients that capture the airfoil shapes
- `generate_gif.py`--Used to generate the run history of the Genetic Algorithm.  Also produces a rainbow plot similar to `compare_airfoils`
- `rainbow_plot_with_comparison.py`--Allows a simple addition of a comparison airfoil to the rainbow plot

The `postprocess.py` script is a legacy script used for cases before Case 73 (ie, the original OSO WT2 airfoils).  This script must be copied into the run directory of interest and then run using:
```python
mpirun -n 8 python -m mpi4py postprocess.py 0 1
```