	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: September 30, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: October 05, 2015
% Last Modified time: 20:59:37


% Purpose: Simple Kriging with Matern covariance function



function [pred_notrend, sig2pred] = spKrigingMatern(z, locs, pred_locs, nu, rho, sig2, sig2eps)
%% Input Arguments:
% z: detailed residuals
% locs: locations of observations
% pred_locs: prediction locations
% nu: smoothness parameter
% sig2: marginal variance
% sig2eps: variance of measurement error

n = length(z);
%distmat = zeros(n,n);
%for i=1:n
%   distmat(:,i) = great_circle_distance(locs, locs(i,:));
%end
distmat = pdist2(locs,locs);

corr_mat = matern(distmat, nu, rho);
covmat = sig2*corr_mat + sig2eps*speye(n);
L = chol(covmat, 'lower');

m = size(pred_locs,1);
%dist0 = zeros(n,m);
%for i=1:m
%   dist0(:,i) = great_circle_distance(locs, pred_locs(i,:));
%end
dist0 = pdist2(locs, pred_locs);

cross_mat = sig2*matern(dist0, nu, rho);
%E = sparse(n,m);
%for i=1:m
%	E(:,i)=pred_locs(i,1)==locs(:,1) & pred_locs(i,2)==locs(:,2);
%end

%cross_mat = cross_mat + sig2eps*E;

pred_notrend = cross_mat.'*(L'\(L\z));



if nargout > 1
LSig = L\cross_mat;
sig2pred = sig2 - sum(LSig.*LSig,1).';

end
