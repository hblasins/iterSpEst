1. Introduction
---------------
This repository contains the Matlab implementation of the ADMM spectral
estimation agorithm described in

H. Blasinski, J. Farrell and B. Wandell; 'Iterative algorithm for spectral
estimation with spatial smoothing,' ICIP 2015, Quebec City

Please cite the above paper if you use this code in your work.


2. Additional dependencies
--------------------------
The ADMM algorithm uses functions from Matlab's image processing toolbox,
such as the image convolution. Some of the test scripts use the cvx 
toolbox, which can be downloaded from: http://cvxr.com

Certain scripts from the ISET scripts directory require the Image Systems
Engineering Toolbox. They are provided for reference, since all their
outputs are provided in the Data directory.

Some scripts take advantage of the parallel processing toolbox, but they
be trivially adapted to single thread execution.


3. Directory structure
----------------------

a. Data
The data directory has example image data used as inputs to the algorithms. 

b. Figures
This directory contains the scripts used to generate the figures from
the paper. You will need to generate the result files first.

c. Results
The Results directory will contain calculation results that will be produced
by the s_analyzeXXX.m scripts from the root directory.

d. ISET scripts
A directory containing ISET camera simulation scripts. If you have an ISET
licence you can run those scripts to re-create some of the data in the Data
directory.

e. Parameters
A directory containing various parameters such as illuminant SPDs, Macbeth chart
reflectances etc.

e. Utilities
A variety of functions used throughout this project. 

4. Scripts
----------
The test_xxx.m scripts are created to compare the results of the ADMM solver
and the general purpose CVX solver.

The s_analyzeXX.m scripts run the ADMM algorithm on the provided data and 
save the results in the Results directory.

The s_createXXX.m scripts generate figures and tables that were included
in the paper. These figures are stored in the Figures directory. 

The s_simulateXXX.m scripts are ISET simulation scripts.
