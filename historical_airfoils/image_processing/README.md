Image Processing
----------------

This folder shows the process of fitting airfoils from png files.  

- The original airfoil files are located in the `png_files` folder
- The most up to date script is the `process_airfoils.py` script
    - Images should not have a vertical line connecting upper and lower trailing edges
    - This was used for the MHKF1 hydrofoils
- The `process_airfoils_with_green_dot.py` script was previously used to process the wind turbine airfoils, where a green dot was manually placed at the desired (0,0)
    - This is now a legacy script
- The `verification_plots` folder stores some helpful debugging plots
- The `fitted_datfiles` folder contains the raw fits to the data
- The `corrected_datfiles` folder contains the final dat files, which have been thickness corrected