Wind Turbine Airfoil Design Tools for Open-Source Offshore (OSO) Airfoils
=========================================================================

This repository contains the design tools used to produce the Open-Source Offshore (OSO) Airfoils.  This project was led by Sandia National Laboratories, in collaboration with Californa State University, Long Beach (CSULB) and the National Renewable Energy Laboratory (NREL).

This repository is intentionally not set up as a python package for maximum flexibility when running on different computing architectures and systems.  Therefore, you will see such files as `kulfan.py` replicated in many places in this repository.  These are all the same file, but the presence of multiple files allows for the folder to be copy-pasted without the need to install a package.

Note that a `git clone` on this repository will ignore most of the data located in the `postprocessing` folder by design to limit the size of the distribution.


Installation
------------

At present, this code has been tested on MacOS and on a Windows machine running WSL (eg, Linux).  Native Windows support is not currently expected nor guaranteed.  Windows users are advised to set up WSL, at which point the code should run with no issues.

After installing dependencies, the airfoil optimization may be run simply by following the directions in the `runfiles` directory.

The following is a quickstart:
```
git clone <appropriate_ssh_or_html>
cd runfiles
mpirun -n 8 python -m mpi4py c64_t21_l15_r122_k16_n400.py
```

On a 2022 M1 MacBook Air, this runs roughly 800 generations in 48 hours.


Dependencies
------------

The use of these tools assumes the following dependencies, all of which should be installable with `pip install <package>` or `conda install <package>` as appropriate.

- `numpy`
- `scipy`
- `pandas`
- `matplotlib`
- `mpi4py`
- `pint`

We also assume that there is an `xfoil` executable located somewhere in your path.  EG: if you type `which xfoil` in a terminal, a pathname should be printed.  

You may choose to compile XFOIL on your own from source (https://web.mit.edu/drela/Public/web/xfoil/), however, we recommend simply obtaining XFOIL through a distribution of Engineering Sketch Pad (https://acdl.mit.edu/ESP/) (readme is here: https://acdl.mit.edu/ESP/ESPreadme.txt)

Citations
---------

If referencing this work, please cite the following paper:

```
@inbook{karcher2025design, 
    author={Karcher, Cody J. and Maniaci, David C. and Kelley, Chris and Hsieh, Alan and deVelder, Nathaniel and Gupta, Anurag}, 
    title={Design of a preliminary family of airfoils for high reynolds number wind turbine applications}, 
    booktitle={AIAA SCITECH 2025 Forum}, 
    publisher={American Institute of Aeronautics and Astronautics}, 
    place={Orlando, Florida}, 
    month={Jan},
    year={2025},
    DOI={10.2514/6.2025-0840},
    ISBN={978-1-62410-723-8}
    }
```

A pdf of this paper is included in the `publications` folder.


License
-------

Use and distribution of this work is subject to the included MIT License.


Copyright
---------

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.


Funding
-------

We are grateful for the funding that made this work possible provided by the U.S. Department of Energy Office of Energy Efficiency and Renewable Energy Wind Energy Technologies Office.