% Author: Pulong Ma <mpulong@gmail.com>
% Date: December 16th, 2014
% Last Modified by: Pulong Ma
% Last Modified Date: May 11th, 2015
% Last Modified time: 21:03:36


% Purpose: 
% compute the posterior mean of eta given the data z 
% under two different models: 
% 1. SRE+homogeneous covariance model
% 2. SRE+CAR 



function etahat = postmean_eta(z, S, K, Sig_xi, car)
%% input: 
% z is detail residuals
% S is basis function matrix
% K is the covariance matrix of eta in small-scale variation
% Sig_xi is the covariance matrix of xi when car=0;
% is the precision matrix of xi when car=1;
% car is 1 when the model of xi is CAR model, 
% and 0 when the model of xi is multiple of identity.

%% output:
% etahat is the posterior mean of eta given the data z


if nargin < 5
	car = 0;
end

switch car
case 0
	r = size(K, 1);
	n = size(S, 1);
	Sig_xiInv = spdiags(diag(Sig_xi).^(-1), 0, n, n);
	temp = K\speye(r) + S'*(Sig_xiInv*S); 
	SigInv2 = temp\(S'*(Sig_xiInv*z));
	KSDInv = K*(S'*(Sig_xiInv*z));
	etahat = KSDInv - K*S'*(Sig_xiInv*(S*SigInv2));
case 1
	Q = Sig_xi;	
	KSQ = K*S'*Q;
	temp = K\speye(size(S,2)) + S'*Q*S;
	etahat = KSQ*z - KSQ*S*(temp\(S'*(Q*z)));
end

