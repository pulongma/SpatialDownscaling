% Author: Pulong Ma <mpulong@gmail.com>
% Date: December 20th, 2014 
% Last Modified by: Pulong Ma
% Last Modified Date: April 12th, 2015
% Last Modified time: 22:58:10


% Purpose: 




%%%%%%%%%%%%%%%   fixed rank kriging   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pred, sig2FRK]=FRK(data, pred_locs, K, sigxi, sig2_eps,...
						 basis_center, threshold)
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


% observed data
lat = data(:,1);
lon = data(:,2);
z = data(:,3);

% coordinates of prediction locations
lat_pred = pred_locs(:,1);
lon_pred = pred_locs(:,2); 

n = size(z,1);  % number of observations
m = length(lon_pred); % number of prediction locations


% BFs evaluated at observed locations
So = localbasis([lat, lon], basis_center, threshold);

r = size(So, 2);
% BFs evaluated at prediction locations
Sp = localbasis([lat_pred, lon_pred], basis_center, threshold);


%%%%%%%%%%  smooth the data  %%%%%%%%%%%%%%%

% diagonal error variance matrix
D=sigxi*speye(n) + sig2_eps*speye(n);
DInv=inv(D);

% calculate r x r part of the inverse of Sigma
% temp=inv(inv(K)+So'*DInv*So);
temp = (K\speye(r)+So'*DInv*So)\speye(r);

% indicator matrix for observed locations
E = sparse(m,n);
for i = 1:n,
	E((lat_pred==lat(i) & lon_pred==lon(i)), i) = 1;
end

% predict smoothed residuals
SigInvZ = DInv*z - DInv*(So*(temp*(So'*(DInv*z))));
pred = Sp*(K'*So'*SigInvZ) + sigxi*E*SigInvZ;



%%%%%%%%  calculate the FRK variance  %%%%%%%%
if nargout > 1
t0=tic;    
% find rows that correspond to observed locations
index = zeros(m,1);
[indrow, indcol] = find(E==1);
for j = 1:length(indrow)
    index(indrow(j), 1) = indcol(j); 
end

% part 1
SpK = Sp*K;
p1 = sum(SpK.*Sp, 2);
%{
p1 = zeros(m,1);
for i = 1:m
    p1(i,1) = SpK(i,:)*Sp(i,:)';
end
%}

% part 3
KSSigInv = K'*So'*DInv-(K'*So'*DInv*So*temp)*So'*DInv;
SpKSSigInvSK = Sp*(KSSigInv*So*K);
SigInv1 = DInv*So*temp;
SigInv2 = So'*DInv;
p3 = zeros(m,1);
for i = 1:m,
    p3(i, 1) = SpKSSigInvSK(i, :)*Sp(i, :)'; 
    ind = index(i, 1);
    if ind > 0,
        ESigInvE = DInv(ind, ind)-SigInv1(ind, :)*SigInv2(:, ind);
        p3(i,1) = p3(i,1) + 2*sigxi*Sp(i, :)*KSSigInv(:, ind) + sigxi^2*ESigInvE; 
    end;
end

%%% put the pieces of variance together
sig2FRK=p1 + sigxi-p3;
t=toc(t0);
%disp(strcat('time to compute kriging variance:', num2str(t/60), 'mins'))

end
