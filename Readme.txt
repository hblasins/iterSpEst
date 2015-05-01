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


3. Directory structure
----------------------
This root directory constains three subdirectories

a. Data
The data directory has example image data and camera model needed to run the
algorithm. 

b. Figures and Results
These sub-directories are initially empty, and are populated by running 
appropriate scripts. The Results directory will contain calculation results and
the Figures directory will store .eps files with figures used in the paper.


4. Scripts
----------
The test_xxx.m scripts are created to compare the results of the ADMM solver
and the general purpose CVX solver.

The s_analyzeXX.m scripts run the ADMM algorithm on the provided data and 
save the results in the Results directory.

The s_createXXX.m scripts generate figures and tables that were included
in the paper. These figures are stored in the Figures directory. 
