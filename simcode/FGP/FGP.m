% Author: Pulong Ma <mpulong@gmail.com>
% Date: August 31, 2015 
% Last Modified by: Pulong Ma 	
% Last Modified Date: February 11, 2016
% Last Modified time: 22:52:48


% Purpose: Spatial prediction in SRE + CAR model for all the 
% locations (including observed and unobserved)

function [pred, sig2SIK] = FGP(data, pred_locs, K, tau2, ...
		gamma, H, A, basis_center, threshold, sig2eps, Delta)
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

% observed data
lat = data(:,1);
lon = data(:,2);
z = data(:,3);


% coordinates of prediction locations
lat_pred = pred_locs(:,1);
lon_pred = pred_locs(:,2); 


n = size(z,1);  % number of observations
m = length(lon_pred); % number of prediction locations
N = size(H, 1);
Veps = sig2eps*speye(n);

if nargin < 11 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
else % weighted CAR
	DeltaInv = spdiags(diag(Delta).^(-1), 0, N, N);
end	

% BFs evaluated at observed locations
So = localbasis([lat, lon], basis_center, threshold);

% BFs evaluated at prediction locations
Sp = localbasis([lat_pred, lon_pred], basis_center, threshold);
r = size(So, 2);

% link matrix
Ao = A{1};
Ap = A{2};

% help terms
VInv = spdiags((diag(Veps)).^(-1), 0, n, n);
AVA = Ao'*VInv*Ao;
AVz = Ao'*VInv*z;
AVS = Ao'*(VInv*So);

%%%%%%%%%%  smooth the data  %%%%%%%%%%%%%%%
Q = (DeltaInv - gamma*H)/tau2;
mid = Q + AVA;
[L1, ~, s1] = lchol(mid);
mid4S(s1, :) = L1'\(L1\AVS(s1, :));
DSo = VInv*So - VInv*(Ao*mid4S); clear mid4S;
SDS = So'*DSo;
temp = K\speye(r) + SDS;
%mid3z(s1, 1) = L1'\(L1\AVz(s1,1));
mid3z(s1, 1) = cs_ltsolve(L1, L1\AVz(s1,1));
Dz = VInv*z - VInv*(Ao*mid3z);
SigInvZ = Dz - DSo*(temp\(So'*Dz));

% predict smoothed residuals
part1 = Sp*(K*(So'*SigInvZ));
[Lq, ~, sq] = lchol(Q);
ASigInvZ = Ao'*SigInvZ;
%QASigInvZ(sq, 1) = Lq'\(Lq\ASigInvZ(sq,1));
QASigInvZ(sq, 1) = cs_ltsolve(Lq, Lq\ASigInvZ(sq,1));

part2 = Ap*QASigInvZ;
pred = part1 + part2;


%%%%%%%%  calculate the SIK variance  %%%%%%%%
if nargout > 1  % upon request
t0 = tic;	
% part1: Sp*K*So'*D*So*K*Sp - SpK*So'*D*So*temp^(-1)*SDS*SpK'
SpK = Sp*K;
p0 = sum(SpK.*Sp, 2);
tempSDS = temp\SDS;
p1 = sum((SpK*SDS*(speye(r)-tempSDS)).*SpK, 2);

% part2: SpK*So'*D*Ao*Q^(-1)*Ap' + Ap*Q^(-1)*Ao'*DSo*SpK'
%  -[SpK*SDS*tempSo'*D*Ao*Q^(-1)*Ap'+ Ap*Q^(-1)*Ao'*DSo*tempSDS*SpK']
ADS = Ao'*DSo; clear DSo;
QInvADS(sq, :) = Lq'\(Lq\ADS(sq, :)); clear ADS;
ApQInvADS = Ap*QInvADS; clear QInvADS;
p2 = 2*sum((ApQInvADS*(speye(r)-tempSDS)).*SpK, 2);
clear SpK;

% part4: Ap*Q^(-1)*Ao'*D*So*temp^(-1)*So'*D*Ao*Q^(-1)*Ap';
p4 = sum((ApQInvADS/temp).*ApQInvADS, 2);
clear ApQInvADS;

% part3: Ap*Q^(-1)*Ao'*D*Ao*Q^(-1)*Ap'
p3 = zeros(m,1);
sig2xi = zeros(m, 1);
for i = 1:m
	QAp(sq, 1) = cs_ltsolve(Lq, Lq\full(Ap(i,sq)'));
	sig2xi(i) = Ap(i,:)*QAp;
	AVAQAp = AVA*QAp;
	mid4QAp(s1, 1) = cs_ltsolve(L1, L1\AVAQAp(s1,1));
	DAQAp = VInv*(Ao*(QAp-mid4QAp));
	p3(i) = (Ao*QAp)'*DAQAp;
end


% putting all pieces together
sig2SIK = p0+sig2xi - (p1+p2+p3-p4);
t=toc(t0);
disp(strcat('time to compute kriging variance:', num2str(t/60), 'mins'))
end



end % end main function


