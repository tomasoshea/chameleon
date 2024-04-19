# solar chameleon production

Code to calculate the solar Primakoff production of light scalars / chameleons.

- *filter.py* reads solar model from AGSS09 files in data folder and outputs easy to use .dat files
- *solarcham.cpp* contains all the code needed to calculate spectra, energy loss, profiles etc in the form of .dat files
- *utils.h* contains constants and that and is used by *solarcham.cpp*
- the python plotting codes then produce various plots based on the data outputted by *solarcham.cpp*
