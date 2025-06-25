Running the Optimization
------------------------

This folder contains all of the relevant files necessary to create your own airfoils or to replicate the work performed in this project.

A typical workflow is as follows:

1. Use the `copy_scripts.py` file to copy the `common_runner.py` script to your desired filenames (see the section below on file keycodes)
2. Run the optimization via the following command: `mpirun -n 188 python -m mpi4py c63_t15_l15_r124_k8_n200.py` or similar
    - Note that fewer processes are recommended for standard laptops or personal computers:  `mpirun -n 8 python -m mpi4py c63_t15_l15_r124_k8_n200.py`
    - A reasonable personal laptop or computer can run a K=8, N=200 case in roughly 24-48 hours and produces a reasonable result


File Keycodes
-------------

Running files must abide by a specific keycode pattern.  Note that keys need not be in a specific order, and this can be easily modified in the run script if desired

- `c##` (required) -- The case number.  Tested as a 2 digit number (eg `01`-`99`), but may work for other numbers as well
- `t##` (required) -- The thickness to chord ratio (commonly known as 'tau') for the airfoil being designed.  At present, must be one of [15, 18, 21, 24, 27, 30, 33, 36]
- `k##` (required) -- The number of design variables (eg Kulfan coefficients), split evenly between top and bottom surfaces.  K=8 is 4 coefficients for the top surface, and 4 for the bottom.  In general, K should be at least 8, but limited gains will be observed at K>16.  This work does not use a K of more than 16 for any of the runs.  Increasing this number too high also allows airfoils to become too tailored to the specific design point and worsens the structural integration challenges.
- `n###` (required) -- The number of airfoils to use in the population.  Recommended value is at least 20 x (K), or 20 x (Number of design variables).  Generally roughly 400 is used for a run of K=16 and 200 is used for K=8.
- `l##` (optional) -- The lift coefficient to use in the design case, divided by 10.  For example, `l10` is a lift coefficient of 1.0, `l17` is 1.7 and so on.  Default is determined by the table below.
- `r(#)##` (optional) -- The lift to drag ratio that must be achieved in the 'rough' analysis case.  May be a two digit number (eg `r98`) or a three digit number (eg `r125`).  Default is determined by the table below.
- `e(#)#` (optional) -- The Reynolds number (in millions) to run the case at.  `e1` would be one million, `e6` six million, `e12` twelve million, `e22` twenty-two million, etc.  Code has only been tested up to 18 million.  Default is determined by the table below.
- `i###` (optional, default=1200) -- The number of generations (or iterations) to run the genetic algorithm for.  Default is 1200.