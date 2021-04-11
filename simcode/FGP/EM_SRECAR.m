	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: August 31, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: July 24, 2017
% Last Modified time: 14:58:24


% Purpose: EM algorithm to estimate parameters in SRE + CAR model
% In this function, only eta is considered to be missing data, and EM 
% algorithm is adopted to estimate parameters K, tau2, gamma, and possible
% sig2eps. For measurement error sig2eps, sometimes sig2eps is known from
% the experiment or data, so there is no need to estimate sig2eps in EM 
% algorithm; while there is a possibility that we do not know sig2eps in 
% advance, or we want to compare the MLE of sig2eps with the true measurement
% error, so we can also estimate sig2eps in EM algorithm along with other 
% parameters.
% This routine does not estimate trend coefficients in EM algorithm!!!

function [K_em, tau2_em, gamma_em, T, epsilon, time] = EM_SRECAR(S, ...
			A, z, H, sig2eps, Delta, maxit, avgtol, opts)
%% input:
% S: is an nxr basis function matrix.
% z: is the detail residuals.
% H: is the proximity matrix for the entire ocean with dimension (m+n)x(m+n)
% sig2eps: is the measurem error. When it is not given, then the function 
% can estiamte this measurement error; otherwise this algorithm will use 
% the true value of sig2eps given by the user.
% Delta: is the (m+n)x(m+n) diagonal matrix
% maxit: is the maximum number of iterations in EM algorithm.
% avgtol: the tolerance error when estimating parameters.
% opts: is an optimzation object specifying the optimzation algorithm.

n = length(z); % number of observed data
r = size(S, 2); % number of basis function
N = size(H,1); % number of area units in discretized domain

% default values for the last two parameters
if nargin<8, avgtol=1e-3; end
if nargin<7, maxit=100; end
if nargin<6 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
	reord = symrcm(H);
	large = eigs(H(reord,reord), 1, 'LA');
	small = eigs(H(reord,reord), 1, 'SA');
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
varest = var(z, 1);
K_old = 0.95*varest*eye(r);
tau2 = 0.05*varest;
t = 1;
done = 0;



update = 0;
lbt = 0.01;
ubt = 10;

lb = [lbt; lbg];
ub = [ubt; ubg];

if nargin < 9
	opts = optimoptions(@fmincon,'TolX', 1e-3, 'Algorithm', ...
		'active-set','Display','off');
end

fun = @neg2Q;

% help term
Veps = sig2eps*speye(N);
AVA = A*Veps*A';
AVAInv = spdiags((diag(AVA)).^(-1), 0, n, n);
A2A = A'*AVAInv*A;
A2z = A'*AVAInv*z;
A2AS = A'*AVAInv*S;

while done == 0
ts(t)=tic;

mid = (DeltaInv-gamma(t)*H)/tau2(t) + A2A;
[L1, ~, s1] = lchol(mid);
mid4S(s1, :) = L1'\(L1\A2AS(s1, :));
DS = AVAInv*S - AVAInv*(A*mid4S);
SDS = S'*DS; clear mid4S DS;
temp = K_old\speye(r) + SDS;
mid3z(s1, 1) = cs_ltsolve(L1, L1\A2z(s1,1)); clear R1;
Dz = AVAInv*z - AVAInv*(A*mid3z);

% update K
muEta = K_old*(S'*Dz) - K_old*(SDS*(temp\(S'*Dz)));
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

% check convergence 
if abs(gamma(t+1)-gamma(t)) + abs(tau2(t+1)-tau2(t)) < 1e-3 
	update = 1;
end

diff = sum(sum((K_new-K_old).^2, 1), 2) + (tau2(t+1)-tau2(t))^2 ...
		+ (gamma(t+1)-gamma(t))^2;
if diff < min(avgtol*r^2,1)
	done = 1;
end

if t > maxit
	done = 1;
	%disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
end
te(t) = toc(ts(t));

K_old = K_new;
t = t+1;
%disp(strcat(' t= ', num2str(t)))
end  % end while


% check positive definiteness of K
E = eig(K_new);
if nnz(E<=0) > 0
	disp('Error: K is NOT Positive definite!')
end

K_em = K_new;
tau2_em = tau2(t);
gamma_em = gamma(t);

if nargout > 3
	T = t;
end

if nargout > 4
	epsilon = diff;
end	

if nargout >5
	time = te;
end

function f = neg2Q(theta)
	Q = (DeltaInv-theta(2)*H)/theta(1);
	Tmid = Q + A2A;
	[Lq, ~, sq] = lchol(Q);
	logQ = 2*sum(log(diag(Lq))); clear Lq sq;
	[Lt, ~, st] = lchol(Tmid);
	logT = 2*sum(log(diag(Lt)));
	part1 = -logQ+logT;
	Tmid4S(st,:) = Lt'\(Lt\A2AS(st, :));	  
	DnewS = AVAInv*S - AVAInv*(A*Tmid4S);
	clear Tmid4S;
	%Tmid3Z(st, 1) = Lt'\(Lt\A2z(st,1)); clear Lt;
	Tmid3Z(st, 1) = cs_ltsolve(Lt, Lt\A2z(st,1)); clear Lt;
	DnewZ = AVAInv*z - AVAInv*(A*Tmid3Z);
	part2 = z'*DnewZ - 2*(z'*DnewS)*muEta;
	SDnewS = S'*DnewS; 
	part3 = sum(sum(SDnewS.*SigEta', 2)) + muEta'*(SDnewS*muEta);
	f = part1 + part2 + part3;

end  % end nested function



end
