	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: May 12th, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: May 12th, 2015
% Last Modified time: 


% Purpose: 


function [xihat, sig2_xi_hat] = posterior_xi(z, S, K, Q, car)
%% input:
% z is detail residuals
% S is basis function matrix
% K is the covariance matrix of eta in small-scale variation
% Q is the covariance matrix of xi in fine-scale variation    
% when car=0; is the precision matrix of xi when car=1;
%% output:
% xihat is the posterior mean of xi given the data z
% Sig_xi_hat is the posterior variance of xi given the data z


if nargin < 5
	car = 1;
end

switch car 
case 0
	var_xi = diag(Q);
	D = spdiags(var_xi.^(-1), 0, size(Q));
	temp = K\speye(size(S,2)) + S'*D*S;
	xihat = z - S*(temp\(S'*(D*z)));
	if nargout > 1
		Sig_xi_hat = S*(temp\S');
	end	
case 1
	temp = K\speye(size(S,2)) + S'*Q*S;
	xihat = z - S*(temp\(S'*(Q*z)));
	if nargout > 1
		n = size(S, 1);
		sig2_xi_hat = zeros(n, 1);
		parpool('local', 4);
    	parfor i = 1:n
    		sig_xi_hat(i) = (S(i, :)/temp)*S(i,:)';
    	end
    	delete(gcp('nocreate'));
	end
end



