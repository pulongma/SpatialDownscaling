

clear; clc; 


addpath('./code')
addpath('./FGP')


f = @(x1, x2) 50*x1.*exp(-x1.^2-x2.^2)

nx = 100;
ny = 100;
s1 = -2;
s2 = 6;
L0 = s2 - s1;
xg = linspace(s1, s2, nx)';
yg = linspace(s1, s2, ny)';
[x, y] = meshgrid(xg, yg);
x = reshape(x, nx*ny, 1);
y = reshape(y, nx*ny, 1);
total_locs = [x, y];

F = f(x,y);
sig2eps = 0.01*var(F);

rng(2345);
Y = F + sqrt(sig2eps)*randn(nx*ny,1);

%%% prediction locations 
n = nx*ny;
m = floor(n*0.1);

flagMR = total_locs(:,1)>=-0.5 & total_locs(:,1)<=0.5 & ...
 total_locs(:,2)>=-1.5 & total_locs(:,2)<=1.5;


total_ID = 1:n;
total_ID = total_ID(:); 
MBD_ID = total_ID(flagMR==1);
RemainID = total_ID(flagMR==0);


rng(230);
indrand = datasample(RemainID, m, 'Replace', false);
ind_locs = setdiff(RemainID, indrand);


MBD = MBD_ID;
MAR = indrand;
predID = total_ID;
pred_locs = total_locs(predID, :);
locs = total_locs(ind_locs, :);

Z = Y(ind_locs);


%%%%%%% construct adaptive basis functions
load ./basisfun_cutoff90_sim new_threshold new_center

basis_center = new_center{12};
radius = new_threshold{12};

S = localbasis(locs, basis_center, radius);

Ao = indicator_mat(locs, total_locs);
dtemp = pdist2(total_locs, total_locs);
[indi, indj] = find(dtemp>0 & dtemp<0.081); % first neighbor hood
H = sparse(indi, indj, ones(length(indi),1), size(total_locs,1), size(total_locs,1));
reord = symrcm(H);
large = eigs(H(reord,reord), 1, 'LA');
small = eigs(H(reord,reord), 1, 'SA');
Hall.H = H;
Hall.large = large;
Hall.small = small;

%%%%%%% EM estimation 
Xo = ones(size(Ao,1), 1);
%[K,tau2,gamma]=EM_SRECAR(S,Ao,Z,H,sig2eps);
tic;
[K,tau2,gamma,beta]=EM_FGP(S,Ao,Z,Xo,Hall,sig2eps);
Q =  (speye(size(H,1)) - gamma*H)/tau2;


%%%%%%% prediction
data = [locs, Z];
pred_locs = total_locs;
Ap = indicator_mat(pred_locs, total_locs);
A = {Ao, Ap};
[predSIK] = FGP(data,pred_locs,K,tau2,gamma,H,A,basis_center,radius,sig2eps);
toc;
Xp = ones(size(pred_locs,1), 1);
predSIK = predSIK + Xp*beta;
indmiss = [MAR;MBD];
mean((F(indmiss)-predSIK(indmiss)).^2)
mean((F(MAR)-predSIK(MAR)).^2)
mean((F(MBD)-predSIK(MBD)).^2)


