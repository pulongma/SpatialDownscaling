

clear; clc; 


addpath('./code')
addpath('./FRK')


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

%flagMR = zeros(n,1);

%flagMR = total_locs(:,1)>s1+0.2*L0 & total_locs(:,1)<s1+0.6*L0 & ...
% total_locs(:,2)>s1+0.2*L0 & total_locs(:,2)<s1+0.5*L0;

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


% define initial basis-function
a=5;
b=5;
cx = linspace(-2, 6, a);
cy = linspace(-2, 6, b);
[cx,cy] = meshgrid(cx, cy);

basis1 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis1);
radius1 = 1.5*min(d)*ones(a*b, 1);


a=7;
b=7;
cx = linspace(-2, 6, a);
cy = linspace(-2, 6, b);
[cx,cy] = meshgrid(cx, cy);

basis2 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis2);
radius2 = 1.5*min(d)*ones(a*b, 1);


a=9;
b=9;
cx = linspace(-2, 6, a);
cy = linspace(-2, 6, b);
[cx,cy] = meshgrid(cx, cy);

basis3 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis2);
radius3 = 1.5*min(d)*ones(a*b, 1);


a=11;
b=11;
cx = linspace(-2, 6, a);
cy = linspace(-2, 6, b);
[cx,cy] = meshgrid(cx, cy);

basis4 = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis4);
radius4 = 1.5*min(d)*ones(a*b, 1);



grid_all1 = [basis1;basis2;basis3];
threshold1 = [radius1;radius2;radius3];
nS1 = ones(length(threshold1),1);


grid_all2 = [basis1;basis2;basis3;basis4];
threshold2 = [radius1;radius2;radius3;radius4];
nS2 = ones(length(threshold2),1);


FLAG1_basis=grid_all1(:,1)>=-0.5 & grid_all1(:,1)<=0.5 & ...
 grid_all1(:,2)>=-1.5 & grid_all1(:,2)<=1.5;


FLAG2_basis=grid_all2(:,1)>=-0.5 & grid_all2(:,1)<=0.5 & ...
 grid_all2(:,2)>=-1.5 & grid_all2(:,2)<=1.5;


S = localbasis(locs, grid_all1(FLAG1_basis==0,:), threshold1(FLAG1_basis==0));

tic;
[K1, sig2_xi1, T] = EM_FRK(S,Z,sig2eps);

pred_locs = total_locs;

[pred,sig2FRK] = FRK([locs, Z], pred_locs, K1, sig2_xi1, sig2eps, grid_all1(FLAG1_basis==0,:), threshold1(FLAG1_basis==0));
toc;

indmiss = [MAR;MBD];
mean((Y(indmiss)-pred(indmiss)).^2)
mean((Y(MAR)-pred(MAR)).^2)
mean((Y(MBD)-pred(MBD)).^2)



S = localbasis(locs, grid_all2(FLAG2_basis==0,:), threshold2(FLAG2_basis==0));

tic;
[K2, sig2_xi2, T] = EM_FRK(S,Z,sig2eps);

pred_locs = total_locs;
[pred,sig2FRK] = FRK([locs, Z], pred_locs, K2, sig2_xi2, sig2eps, grid_all2(FLAG2_basis==0,:), threshold2(FLAG2_basis==0,:));

toc;


indmiss = [MAR;MBD];
mean((Y(indmiss)-pred(indmiss)).^2)
mean((Y(MAR)-pred(MAR)).^2)
mean((Y(MBD)-pred(MBD)).^2)





