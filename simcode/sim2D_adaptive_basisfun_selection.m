
clear; clc; 

addpath('./autobasis')
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


%%%%% figure 
Z_old =NaN(nx*ny, 1);
Z_old(ind_locs) = Z;

figure(6); clf;
clims = [ceil(min(Z_old)), floor(max(Z_old))];
fimage=imagesc(xg, yg, reshape(Z_old,nx,ny), clims);
rectangle('Position', [-0.5,-1.5,1,3], 'LineWidth', 1);
set(gca,'YDir','normal')
%set(gcf,'AlphaData', ~isnan(Z_old))
title('Simulated Data')
xlim([s1, s2]);
ylim([s1, s2]);

nanpa=[1 1 1;parula];
colormap(gcf,nanpa)
%dmap=(max(clims)-min(clims))/size(parula,1);
%caxis([ceil(min(clims)-dmap), ceil(max(clims))])
hcor=colorbar;
ylim(hcor, clims)


figure(7); clf;
clims = [ceil(min(Z_old)), floor(max(Z_old))];
imagescwithnan(xg, yg, reshape(Z_old,nx,ny), parula, [1 1 1]);
rectangle('Position', [-0.5,-1.5,1,3], 'LineWidth', 1);
set(gca,'YDir','normal')
%set(gca,'AlphaData', ~isnan(Z_old))
title('Simulated Data')
xlim([s1, s2]);
ylim([s1, s2]);
colorbar;


% define initial basis-function
a=5;
b=5;
cx = linspace(-2, 6, a);
cy = linspace(-2, 6, b);
[cx,cy] = meshgrid(cx, cy);

basis = [reshape(cx, a*b, 1), reshape(cy, a*b, 1)];
d = pdist(basis);
radius = 1.5*min(d)*ones(a*b, 1);
nS = ones(a*b, 1);

% calculate residuals with initial basis functions

grid_all = basis;
threshold = radius;
nS = ones(length(threshold), 1);
S = localbasis(locs, grid_all, threshold);
[Ssub, indsub] = licol(S,1);
basis = grid_all(indsub, :);
radius = threshold(indsub);


grid_all = basis;
threshold = radius;
nS = ones(length(threshold), 1);
S = localbasis(locs, grid_all, threshold);
[K, sig2_xi, T] = EM_FRK(S,Z,sig2eps);


[res_old] = BPC(locs, Z, grid_all, threshold, sig2eps);



data{1} = locs(:, 1);
data{2} = locs(:, 2);
data{3} = res_old;

n = size(locs, 1);

ind = cell(n, 1);
tic;
%parpool('local', 2);
for i = 1:n
	dtemp = great_circle_distance(locs, locs(i, :));
	dtemp(locs(:,1)==locs(i,1) & locs(:,2)==locs(i,2)) = 0;
	ind{i} = find(dtemp <= 0.5);
end
%delete(gcp('nocreate'));
toc;
data{4} = ind; 


criteria = 5;

i=1;
new_center{i} = basis;
new_threshold{i} = radius;
cutoff(i) = 100;
r = size(new_center{i}, 1);
tol=100;
res_new{i} = res_old;
ELMSE{i} = zeros(n,1);
done=0;

while(done==0)
	tstart=tic;
	display(strcat('iteration i=', num2str(i)));
	[new_center{i+1}, new_threshold{i+1}, cutoff(i+1), res_new{i+1}, stop_update, ELMSE{i}] = basisfun(data, Z, ...
		new_center{i}, new_threshold{i}, criteria, sig2eps);	
	timing(i)=toc(tstart);
	tol = abs(cutoff(i)-cutoff(i+1));
	r = size(new_center{i},1);
	i = i+1;
	data{3} = res_new{i};
	if stop_update==1
		break;
	end
	if tol<0.01
		done=1;
	else
		if r>250 
			done=1;
		else
			done=0;
		end
	end

end
toc;


res_old = data{3};
ind = data{4};
sd_res = zeros(n, 1);
mu_res = zeros(n, 1);
rule = zeros(n, 1);
rule_old = zeros(n, 1);
rule_old = zeros(n, 1);
%parpool('local', 2);
for k = 1:n
	res_temp = res_old(ind{k});
	sd_res(k) = std(res_temp);
	mu_res(k) = mean(res_temp);
	rule_old(k) = (mu_res(k))^2 + sd_res(k)^2;
end
ELMSE{i} = rule_old;


save ./basisfun_cutoff90_sim res_new new_threshold new_center res_new cutoff timing ELMSE



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% figures to diagnose critera measure
iter = 1:12;
figure(1); clf;
fh(1)=subplot(1,2,1);   
boxplot(cell2mat(ELMSE(iter)), 'symbol', 'r.');
xlim([0, max(iter)])
xlabel('iteration')
ylabel('ELMSE')
% creat new pair of axes inside current figure
fh(2)=axes('position', [.25, .565, .19, .29])
box on; 
indexOfInterest = iter>=6 & iter<=max(iter);
labels = strtrim(cellstr(num2str(iter(indexOfInterest)'))')
boxplot(cell2mat(ELMSE(iter(indexOfInterest))), 'symbol', 'r.', ...
	'Labels',labels)
axis tight
ylim([0, 5])

fh(3)=subplot(1,2,2);
plot(iter, cutoff(iter), 'k--^','MarkerSize', 5, 'LineWidth', 1)
xlim([0, max(iter)])
xlabel('iteration')
ylabel('cutoff')

% creat new pair of axes inside current figure
fh(4)=axes('position', [.69, .575, .19, .28])
box on; 
indexOfInterest = iter>=6 & iter<=max(iter);
plot(iter(indexOfInterest), cutoff(indexOfInterest), ...
	'k--^','MarkerSize', 4, 'LineWidth', 1);
axis tight

xlim([6, max(iter)])
ylim([0, 2])

set(findobj('type','axes'),'fontsize',12)
yticks(fh(2),[0, 2, 4])
xticks(fh(4), [6:max(iter)])
xticks(fh(3), iter)

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])

print(gcf, './basis_measures_sim', '-depsc')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_all = new_center{12};
threshold = new_threshold{12};
nS = ones(length(threshold), 1);
S = bisquare(locs, grid_all, threshold);
tic;
%%%%%% FRK with adaptive basis functions
[K, sig2_xi] = EM_FRK(S,Z,sig2eps);

pred_locs = total_locs;
[pred] = FRK([locs, Z], pred_locs, K, sig2_xi, sig2eps, grid_all, threshold);
toc;

indmiss = [MAR;MBD];
mean((F(indmiss)-pred(indmiss)).^2)
mean((F(MAR)-pred(MAR)).^2)
mean((F(MBD)-pred(MBD)).^2)

%%%%% FGP with adaptive basis functions
addpath('./FGP')
Ao = indicator_mat(locs, total_locs);
dtemp = pdist2(total_locs, total_locs);
[indi, indj] = find(dtemp>0 & dtemp<0.081); % first neighborhood structure
H = sparse(indi, indj, ones(length(indi),1), size(total_locs,1), size(total_locs,1));
reord = symrcm(H);
large = eigs(H(reord,reord), 1, 'LA');
small = eigs(H(reord,reord), 1, 'SA');
Hall.H = H;
Hall.large = large;
Hall.small = small;

Xo = ones(size(Ao,1), 1);

tic;
%%%%%%% EM estimation 
%[K,tau2,gamma]=EM_SRECAR(S,Ao,Z,H,sig2eps);
[K,tau2,gamma,beta]=EM_FGP(S,Ao,Z,Xo,Hall,sig2eps);
Q =  (speye(size(H,1)) - gamma*H)/tau2;


%%%%%%% prediction
data = [locs, Z];
pred_locs = total_locs;
Ap = indicator_mat(pred_locs, total_locs);
A = {Ao, Ap};
[predSIK] = FGP(data,pred_locs,K,tau2,gamma,H,A,grid_all,threshold,sig2eps);
toc;
Xp = ones(size(pred_locs,1), 1);
predSIK = predSIK + Xp*beta;

indmiss = [MAR;MBD];
mean((F(indmiss)-predSIK(indmiss)).^2)
mean((F(MAR)-predSIK(MAR)).^2)
mean((F(MBD)-predSIK(MBD)).^2)



