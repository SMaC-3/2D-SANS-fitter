# 2D SANS fit optimiser
Python code to fit 2D, anisotropic small-angle neutron scattering (SANS) data based on sasmodels functionality. Model fits are optimised using the Levenberg-Marquardt algorithm from scipy. Requires the open-source sasmodels package: https://github.com/SasView/sasmodels

## Features

- Fit 1D or 2D SANS data using scattering models available in sasmodels 
- Perform sector and annular extractions 
- Plot experimental and simulated data 

## How to use

2D, anisotropic data fitting is performed using the script rheoSANS_fitOpt.py

This script imports rheoSANS_fitOpt_Functions.py to make use of its custom class (sans2d) and additional functions. 

The following setup steps are required:

- A user must create a csv file with column headings 'index', 'filename', 'sample', and 'shear'. Modifications to the csv and data path will be necessary under rheoSANS_fitOpt_Functions.py > class sans2d > method getData
- The selected sasmodels model can be changed under rheoSANS_fitOpt_Functions.py > class sans2d > method makeCalc
- The code was initially developed for a shear-banding fluid. This refers to a macro scale fracture of a sample during flow, resulting in two regions of with different levels of anisotropy.  A shear band parameter is included to account for this, but should be set to 1 for non-shear banded samples.  
- Fixed and initial guess parameters are supplied through two dictionaries. Only pars2sim is required for non-shear banded samples. For shear-banding samples, the low shear band parameters are fixed and provided in par2static. 
- Set indexSelected to import the desired data according to the reference csv file, and modify the error check parameters
- Use the dictionary fitChoose to select the fitting parameters 
- Run the code to optimise the model. 




