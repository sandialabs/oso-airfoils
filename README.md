Wind Turbine Airfoil Design Tools for Open-Source Offshore (OSO) Airfoils
=========================================================================

This repository contains the design tools used to produce the Open-Source Offshore (OSO) Airfoils.  This project was led by Sandia National Laboratories, in collaboration with Californa State University, Long Beach (CSULB) and the National Renewable Energy Laboratory (NREL).

This repository is intentionally not set up as a python package for maximum flexibility when running on different computing architectures and systems.  Therefore, you will see such files as `kulfan.py` replicated in many places in this repository.  These are all the same file, but the presence of multiple files allows for the folder to be copy-pasted without the need to install a package.


Installation
------------

We strongly recommend a sparse checkout to minimize the required hard drive space, and is achieved via the following:
```
cd <desired_parent_directory>
mkdir oso-airfoils
cd oso-airfoils
git init
git remote add -f origin git@github.com:sandialabs/oso-airfoils.git #use html if appropriate
git config core.sparseCheckout true
echo "historical_airfoils/" >> .git/info/sparse-checkout
echo "publications/" >> .git/info/sparse-checkout
echo "released_designs/" >> .git/info/sparse-checkout
echo "runfiles/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
echo "LICENSE" >> .git/info/sparse-checkout
echo ".gitignore" >> .git/info/sparse-checkout
echo "postprocessing/cached_data" >> .git/info/sparse-checkout
echo "postprocessing/README.md" >> .git/info/sparse-checkout
echo "postprocessing/kulfan.py" >> .git/info/sparse-checkout
echo "postprocessing/postprocess.py" >> .git/info/sparse-checkout
echo "postprocessing/cases/caselog.txt" >> .git/info/sparse-checkout
echo "postprocessing/cases/active/" >> .git/info/sparse-checkout
git pull origin main
git branch --set-upstream-to=origin/main
```

This will take some time to index all of the files (particularly on the step that adds the origin), but will not clone any of the data files onto your hard drive.

A normal clone is still possible:
```
git clone <html_or_ssh_link>
```
but will take up significant hard drive space.

At present, this code has been tested on MacOS and on a Windows machine running WSL (eg, Linux).  Native Windows support is not currently expected nor guaranteed.  Windows users are advised to set up WSL, at which point the code should run with no issues.

After installing dependencies, the airfoil optimization may be run simply by following the directions in the `runfiles` directory.

The following is a quickstart:
```
cd runfiles
mpirun -n 8 python -m mpi4py c65_t21_l15_r122_k16_n400.py
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
    title={Design of a Preliminary Family of Airfoils for High Reynolds Number Wind Turbine Applications}, 
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