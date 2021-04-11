	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: August 31, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: July 24, 2017
% Last Modified time: 14:58:24


% Purpose: EM algorithm to estimate parameters in SRE + CAR model
% In this function, only eta is considered to be missing data, and EM 
% algorithm is adopted to estimate parameters beta, K, tau2, gamma, and possible
% sig2eps. For measurement error sig2eps, sometimes sig2eps is known from
% the experiment or data, so there is no need to estimate sig2eps in EM 
% algorithm; while there is a possibility that we do not know sig2eps in 
% advance, or we want to compare the MLE of sig2eps with the true measurement
% error, so we can also estimate sig2eps in EM algorithm along with other 
% parameters.
% This routine does estimate trend coefficients in EM algorithm!!!

function [K_em, tau2_em, gamma_em, beta_em] = EM_SIK(S, ...
			A, Z, X, Hall, sig2eps, Delta, maxit, avgtol, opts)
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



H = Hall.H;
large = Hall.large;
small = Hall.small;


n = length(Z); % number of observed data
r = size(S,2); % number of basis function
N = size(H,1); % number of area units in discretized domain

% default values for the last two parameters
if nargin<9, avgtol=1e-3; end
if nargin<8, maxit=100; end
if nargin<7 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
	lbg = 1/small + 1e-3*abs(1/small);
	ubg = 1/large - 1e-3*abs(1/large);
	gamma = (lbg + ubg)/2;
else % weighted CAR
	DeltaInv = spdiags(diag(Delta).^(-1), 0, N, N);
	lbg = -1+1e-3;
	ubg = 1-1e-3;
	gamma = 0.16;
end	

% initial values
varest = var(Z, 1);
K_old = 0.95*varest*eye(r);
tau2 = 0.05*varest;
t = 1;
done = 0;
neg2loglikelihood = NaN;

beta_old = regress(Z,X);


update = 0;
lbt = 0.0001;
ubt = 10;

lb = [lbt; lbg];
ub = [ubt; ubg];

if nargin < 10
	opts = optimoptions(@fmincon,'TolX', 1e-3, 'Algorithm', ...
		'active-set','Display','off');
end

fun = @neg2Q;

% help term
Veps = sig2eps*speye(n);
VInv = spdiags((diag(Veps)).^(-1), 0, n, n);
AVA = A'*VInv*A;
AVS = A'*(VInv*S);

VInvX = VInv*X;
AVInvX = A'*VInvX;

while done == 0
ts(t)=tic;

z = Z - X*beta_old;
AVz = A'*(VInv*z);

mid = (DeltaInv-gamma(t)*H)/tau2(t) + AVA;
[L1, ~, s1] = lchol(mid);
mid4S(s1, :) = L1'\(L1\AVS(s1, :));
DS = VInv*S - VInv*(A*mid4S);
SDS = S'*DS; clear mid4S DS;
temp = K_old\speye(r) + SDS;
mid3z(s1, 1) = cs_ltsolve(L1, L1\AVz(s1,1)); 
Dz = VInv*z - VInv*(A*mid3z);

% update beta
muEta = K_old*(S'*Dz) - K_old*(SDS*(temp\(S'*Dz)));
CX(s1, :) = L1'\(L1\AVInvX(s1, :)); clear L1;
DX = VInvX - VInv*(A*CX);
beta_new = (X'*DX)\(DX'*(Z-S*muEta));

% update K
SigEta = K_old - K_old*SDS*K_old' + K_old*SDS*(temp\SDS)*K_old';
K_new = SigEta + muEta*muEta';

% find numerical solution for tau2, gamma
if update == 0
	theta0 = [tau2(t); gamma(t)];
	theta = fmincon(fun, theta0, [], [], [], [], lb, ub, [], opts);
	tau2(t+1) = theta(1);
	gamma(t+1) = theta(2);
elseif update == 1
	tau2(t+1) = tau2(t);
	gamma(t+1) = gamma(t);
end

%{
% calculate negative twice loglikelihood
neg2loglikelihood(t+1) = neg2loglik(K_new,tau2(t+1),gamma(t+1),beta_new);

diff(t) = neg2loglikelihood(t) - neg2loglikelihood(t+1);

if t>=maxit
	done=1;
	disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
else 
	if t>3 & abs(diff(t))<1e-3;
		done=1;
	end
end
%}

% check convergence 

if abs(gamma(t+1)-gamma(t)) + abs(tau2(t+1)-tau2(t)) < 1e-3 
	update = 1;
end

diff = sum(sum((K_new-K_old).^2, 1), 2) + (tau2(t+1)-tau2(t))^2 ...
		+ (gamma(t+1)-gamma(t))^2 + sum((beta_new-beta_old).^2);
if diff < min(avgtol*r^2,1)
	done = 1;
end


if t > maxit
	done = 1;
	%disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
end


te(t) = toc(ts(t));

beta_old = beta_new;
K_old = K_new;
t = t+1;
%disp(strcat(' t= ', num2str(t)));
end  % end while


% check positive definiteness of K
E = eig(K_new);
if nnz(E<=0) > 0
	disp('Error: K is NOT Positive definite!')
end

K_em = K_new;
tau2_em = tau2(t);
gamma_em = gamma(t);

beta_em = beta_new;


function f = neg2Q(theta)
	Q = (DeltaInv-theta(2)*H)/theta(1);
	Tmid = Q + AVA;
	[Lq, ~, sq] = lchol(Q);
	logQ = 2*sum(log(diag(Lq))); clear Lq sq;
	[Lt, ~, st] = lchol(Tmid);
	logT = 2*sum(log(diag(Lt)));
	part1 = -logQ+logT;
	Tmid4S(st,:) = Lt'\(Lt\AVS(st, :));	  
	DnewS = VInv*S - VInv*(A*Tmid4S);
	clear Tmid4S;
	%Tmid3Z(st, 1) = Lt'\(Lt\AVz(st,1)); clear Lt;
	Tmid3Z(st, 1) = cs_ltsolve(Lt, Lt\AVz(st,1)); clear Lt;
	DnewZ = VInv*z - VInv*(A*Tmid3Z);
	part2 = z'*DnewZ - 2*(z'*DnewS)*muEta;
	SDnewS = S'*DnewS; 
	part3 = sum(sum(SDnewS.*SigEta', 2)) + muEta'*(SDnewS*muEta);
	f = part1 + part2 + part3;

end  % end nested function



function f = neg2loglik(K,tau2,gamma,beta)
	z = Z-X*beta;
	AVz = A'*(VInv*z);
	Q = (DeltaInv-tau2*H)/gamma;
	[Lq,~,sq] = lchol(Q);
	[Lc,~,sc] = lchol(Q+AVA);
	part1 = 2*sum(log(diag(Lc))) - 2*sum(log(diag(Lq)));
	clear Lq;

	CInvz(sc,1) = cs_ltsolve(Lc,Lc\AVz(sc,1));
	Dz = VInv*z - VInv*(A*CInvz);

	CInvS(sc,:) = Lc'\(Lc\AVS(sc,:)); clear Lc;
	DS = VInv*S - VInv*(A*CInvS); clear CInvS;
	SDS = S'*DS; clear DS;
	part2 = sum(log(eig(speye(r)+SDS*K)));
	logdetSig = part1 + part2;

	f = logdetSig + z'*Dz - (Dz'*S)*((K\speye(r)+SDS)\(S'*Dz));

end




end
