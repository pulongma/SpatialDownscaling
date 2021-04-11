		
% Author: Pulong Ma <mpulong@gmail.com>
% Date: December 22th, 2014
% Last Modified by: Pulong Ma
% Last Modified date: May 11th, 2015
% Last Modified time: 21:14:24


% Purpose: compute residuals, negative loglikelihood, AIC, BIC in spatial mixed
% effects model





function [res, f, AIC, BIC] = BPC(locs_PCTM, Y, ...
		      grid_all, threshold, sig2eps)

if nargin < 6
	sig2eps=0.25;
end

n = length(Y);


lat = locs_PCTM(:, 1);
lon = locs_PCTM(:, 2);

% create the matrix of predictors (intercept and longitude)
%Xo = [ones(n, 1)];   % constant trend
%Xo = [ones(n, 1), lat];	 
% ordinary least squares estimation of the trend coefficients
%betahat = regress(Y, Xo);
% detrending
%Y_tilde = Y - Xo*betahat;
Y_tilde = Y;


%%%%%%% estimate covariance matrix K of eta %%%%%%%%%%%%%


% SRE + homogeneous covariance model
	% wendland basis function matrix
	S = localbasis(locs_PCTM, grid_all, threshold);
	[K, sig2_xi, T] = EM_FRK(S,Y_tilde,sig2eps);
	Sig_xi = sig2_xi*speye(n);

etahat = postmean_eta(Y_tilde, S, K, Sig_xi);



res = Y_tilde - S*etahat;

if nargout > 1

% Sigma = SKS' + sig2*I
n = size(S, 1);
r = size(S, 2);
V = speye(n);
diagV = diag(V);
DInv = spdiags((sig2_xi*diagV).^(-1), 0, n, n);
temp = K\speye(r) + S'*DInv*S;
%SigInv = DInv - DInv*S*(temp\S')*DInv;
SigInvY = DInv*Y_tilde - DInv*S*(temp\(S'*DInv*Y_tilde));
% matrix determinant lemma
R1 = chol(temp);
R2 = chol(K);
f = 2*sum(log(diag(R1))) + 2*sum(log(diag(R2))) ...
	+ n*log(sig2_xi) + Y_tilde'*SigInvY + n*log(2*pi);
AIC = f + 2*((r^2+r)/2+1);
BIC = f + ((r^2+r)/2+1)*log(n);
%L = exp(f/(-2));

end
