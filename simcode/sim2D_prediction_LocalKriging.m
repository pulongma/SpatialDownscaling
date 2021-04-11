

clear; clc; 


addpath('./code')
addpath('./MVcode')



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
locs_pred = total_locs(predID, :);
locs = total_locs(ind_locs, :);

Z = Y(ind_locs);

%%% prediction with MV approach



tic;
%% prediction at MAR 
m = size(locs_pred(MAR,:),1);
mvsize = 2*(xg(2)-xg(1));

mvk_MAR = zeros(m,1);
sig2mvk = zeros(m,1);
mvrho = zeros(m,1);
mvsig2 = zeros(m,1);

for i=1:m
	[mvk_MAR(i),sig2mvk(i),mvrho(i),mvsig2(i)]=wrapKrige(i,locs_pred(MAR,:),locs,Z,sig2eps,mvsize);
end

MSPE_MAR=mean((F(MAR)-mvk_MAR).^2);

%% prediction at MBD 
m = size(locs_pred(MBD,:),1);
mvsize = 8*(xg(2)-xg(1));

mvk_MBD = zeros(m,1);
sig2mvk = zeros(m,1);
mvrho = zeros(m,1);
mvsig2 = zeros(m,1);

for i=1:m
	[mvk_MBD(i),sig2mvk(i),mvrho(i),mvsig2(i)]=wrapKrige(i,locs_pred(MBD,:),locs,Z,sig2eps,mvsize);
end

MSPE_MBD=mean((F(MBD)-mvk_MBD).^2);


toc;

MSPE_MAR
MSPE_MBD 

res1 = F(MAR)-mvk_MAR;
res2 = F(MBD)-mvk_MBD;

res = [res1;res2];
mean(res.^2)
