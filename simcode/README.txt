####################################################################
This directory contains MATLAB codes to implment basis function selection algorithm as well as make comparison among Local Kriging, Fixed Rank Kriging, and Fused Gaussian Process. For details, see the manuscript titled "Spatial Statistical Downscaling for Constructing High-Resolution Nature Runs in Global Observing System Simulation Experiments"

#######################################################################
The codes for implementing the adaptive basis function selection are in the directory ./autobasis.

The main routine is in the file basisfun.m, which implements the 
adaptive basis function selection algorithm.

%%%% Input arguments:
% 1. data: a cell array containing 
%		   lattitude, 
%		   longitude, 
%          resudals for initial basis function, 
%          index of neighbors for each lattitude/longtitude location
% 2. Y: can be original data or detrended data for corrsponding lat/lon coordinates
% 3. grid_all: a matrix of the centers of basis function in lat/lon form
% 4. threshold: radius of each basis function
% 5. criteria: the cutoff to select basis centers, see below in details
% 6. sig2eps: measurement-error variance 

%%%% Output Arguments:
% 1. new_center: a matrix of new basis centers 
% 2. new_threshold: a vector of new radii for selected basis functions
% 3. cutoff: cutoff of ELMSE
% 4. res_new: a vector of residuals with new basis functions
% 5. stop_update: 1 indicates that no new basis functions are further required.
%                 0 indicates that new basis functions are required.
% 6.elMSE: a vector of empirical-local-mean-squared errors (ELMSE).



######################################################################
The codes for implementing Fixed Rank Kriging are in the directory ./FRK

The main routines are in the files FRK.m and EM_FRK.m

In the file FRK.m, the following inputs and outputs are used.
%%%% Input Arguments:
% 1. data: a matrix with columns containing 
%          lattitude, 
%          longitude, 
%          detrended data 
% 2. pred_locs: a matrix of prediction locations in lat-long form
% 3. K: the r-by-r matrix in the low-rank component.
% 4: sigxi: variance of fine-scale variance term.
% 5. sig2eps: the measurem-error variance. 
% 8. basis_center: a matrix of basis centers in lat-long form.
% 9. threshold: a vector of radii/bandwidths corresponding to basis centers.

%%%% Output Arguments: 
% 1. pred: a vector of prediction mean 
% 2. sig2SIK: a vector of prediction variance


######################################################################
The codes for implementing Fused Gaussian Process are in the directory ./FGP.

The main routines are in the files FGP.m and EM_FGP.m

In the file FGP.m, the following inputs and outputs are used.
%%%% Input Arguments:
% 1. data: a matrix with columns containing 
%		   lattitude, 
%		   longitude, 
%          detrended data 
% 2. pred_locs: a matrix of prediction locations in lat-long form
% 3. K: the r-by-r matrix in the low-rank component.
% 4: tau2: conditional marginal variance in the CAR model.
% 5. gamma: spatial dependence parameter in the CAR model.
% 6. H: an N-by-N proximity matrix
% 7. A: an n-by N basis matrix associated with the CAR model. 
% 8. basis_center: a matrix of basis centers in lat-long form.
% 9. threshold: a vector of radii/bandwidths corresponding to basis centers.
% 10. sig2eps: is the measurem-error variance. 
% 11. Delta: an N-dimensional vector specifying whether homogeneous CAR
% or weighted CAR

%%%% Output Arguments: 
% 1. pred: a vector of prediction mean 
% 2. sig2SIK: a vector of prediction variance


In the file EM_FGP.m, the following inputs and outputs are used.
%%%% Input Arguments:
% 1. S: is an nxr basis function matrix.
% 2. Z: is the data vector.
% 3. X: is the covariates matrix.
% 4. Hall: a structure array containing the proximity matrix and largest and 
% 5. smallest eigenvalues
% 6. sig2eps: is the measurement-error variance. When it is not given, then the function 
% can estiamte this measurement error; otherwise this algorithm will use 
% the true value of sig2eps given by the user.
% 7. Delta: is the NxN diagonal matrix
% 8. maxit: is the maximum number of iterations in EM algorithm.
% 9. avgtol: the tolerance error when estimating parameters.
% 10. opts: is an optimzation object specifying the optimzation algorithm.

%%%% Output Arguments:
% 1. K_em: estimated low-rank covariance matrix K.
% 2. tau2_em: estimated conditional margional variance tau2 in the CAR model.
% 3. gamma_em: estimated spatial dependence parameter gamma in the CAR model.
% 4. beta_em: estimated regression coefficients.


######################################################################
The codes for implmenting Local Kriging are in the directory ./MVcode.


######################################################################
The simulation example to selection basis function adaptively is given in the file sim2D_adaptive_basisfun_selection.m


######################################################################
The prediction based on Local Kriging is in the file sim2D_prediction_LocalKriging.m


######################################################################
The prediction based on equally-spaced basis functions in FRK is in the file sim2D_ES_basis_prediction_FRK.m


######################################################################
The prediction based on adaptive basis functions in FRK is in the file sim2D_adaptive_basisfun_selection.m


######################################################################
The prediction based on adaptive basis functions in FGP is in the file sim2D_adaptive_basisfun_prediction_FGP.m



Copyright (c) [2018] by Pulong Ma <mpulong@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
