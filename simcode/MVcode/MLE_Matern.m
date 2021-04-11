	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: September 28, 2015  
% Last Modified by: Pulong Ma
% Last Modified Date: October 05, 2015
% Last Modified time: 17:52:51


% Purpose: Estimate range parameter rho and marginal variance sig2 in 
% Matern correlation function with maximum likelihood method



function [rho, sig2] = MLE_Matern(z, locs, nu, sig2eps)

%% Input arguments:
% z: detailed residuals
% locs: locations of observations
% nu: smoothness parameter
% sig2eps: variance of measurement error

n = length(z);
%distmat = zeros(n,n);
%for i=1:n
%   distmat(:,i) = great_circle_distance(locs, locs(i,:));
%end
distmat = pdist2(locs,locs);

% initial values
rho_old = 1;
sig2_old = 1;

lb = [0.01; 0.001];
ub = [64/3; 10];



opts = optimoptions(@fmincon,'TolX', 1e-3, 'Algorithm', ...
		'active-set','Display','off');

fun = @neg2Q;

theta0 = [rho_old; sig2_old];
[theta] = fmincon(fun, theta0, [], [], [], [], lb, ub, [], opts);


rho = theta(1);
sig2 = theta(2);


function negloglike = neg2Q(theta);
Corr_mat = matern(distmat, nu, theta(1));
L = chol(theta(2)*Corr_mat + sig2eps*speye(n), 'lower');
negloglike = n*log(2*pi) + 2*sum(log(diag(L))) + z'*(L'\(L\z));
end 


end
