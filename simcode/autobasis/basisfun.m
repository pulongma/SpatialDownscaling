		
% Author: Pulong Ma <mpulong@gmail.com>
% Date: April 15th, 2015
% Last Modified by: Pulong Ma
% Last Modified date: February 10, 2016
% Last Modified time: 22:28:58


% Purpose:  automatically select basis function, i.e., location and radius 
% of each basis function in SRE model

		


function [new_center,new_threshold,cutoff,res_new,stop_update,elMSE] = ...
		basisfun(data, Y, grid_all, threshold, criteria, sig2eps);

%%%% Input arguments:
% 1. data: a cell array containing 
%		   lattitude, 
%		   longitude, 
%          resudals for initial basis function, 
%          index of neighbors for each lattitude/longtitude location
% 2. Y: can be original data or detrended data for corrsponding lat/lon coordinates
% 3. grid_all: a matrix of the centers of basis function in lat/lon form
% 4. threshold: radius of each basis function
% 5. criteria: the cutoff to select basis centers, see below in details
% 6. sig2eps: measurement-error variance 

%%%% Output Arguments:
% 1. new_center: a matrix of new basis centers 
% 2. new_threshold: a vector of new radii for selected basis functions
% 3. cutoff: cutoff of ELMSE
% 4. res_new: a vector of residuals with new basis functions
% 5. stop_update: 1 indicates that no new basis functions are further required.
%                 0 indicates that new basis functions are required.
% 6.elMSE: a vector of empirical local mean-squared errors (ELMSE).

if nargin < 5
	criteria = 2;
end

if nargin < 6
	sig2eps=0.25;
end

locs_PCTM = [data{1}, data{2}];
res_old = data{3};
ind = data{4};

n = size(locs_PCTM, 1);
%ind = cell(n, 1);
sd_res = zeros(n, 1);
mu_res = zeros(n, 1);
rule = zeros(n, 1);
rule_old = zeros(n, 1);
%parpool('local', 2);
for i = 1:n
	res_temp = res_old(ind{i});
	sd_res(i) = std(res_temp);
	mu_res(i) = mean(res_temp);
	rule_old(i) = mu_res(i).^2 + sd_res(i).^2;	
end
%delete(gcp('nocreate'));

switch(criteria)
case 1
	cutoff = quantile(rule_old, 0.7);   
case 2 
	cutoff = quantile(rule_old, 0.75);
case 3
	cutoff = quantile(rule_old, 0.8);
case 4
	cutoff = quantile(rule_old, 0.85);
case 5
	cutoff = quantile(rule_old, 0.9);	
case 6
	cutoff = quantile(rule_old, 0.95);
end

FLAG = find(rule_old>cutoff);
rule_new = rule_old(FLAG);
locs = locs_PCTM(FLAG, :);
res = res_old(FLAG);

disp('Begin updating adaptive basis functions')
r = median(threshold)/3;
h = 0:0.09:(2*r);
delta = 0.1; 

[extrema, eind] = max(rule_new);
center(1, :) = locs(eind, :);

semivari(1, :) = semivariogram(data, center(1, :), r, h, delta);

% delete the center and find the second largest one 
% if it is almost isolated from other data
disp('Delete centers')
notnum = find(isnan(semivari(1, :)));
while (length(semivari(1,:)) - nnz(notnum) < 3)
	rule_new(eind) = eps;
	[extrema, eind] = max(rule_new);
	center(1, :) = locs(eind, :);
	semivari(1, :) = semivariogram(data, center(1, :), r, h, delta);
	notnum = find(isnan(semivari(1, :)));
end



L(1) = LV(2*semivari(1, :), h);
radius(1) = 3*L(1);
separation(1) = 2/3*radius(1);
i = 1;

while (~isempty(locs))									 % separation(i)
	indi = find(great_circle_distance(locs, center(i,:))> separation(i));
	if isempty(indi)
		break
	end
	locs = locs(indi, :);
	rule_new = rule_new(indi);
	[extrema, eind] = max(rule_new);
	center(i+1, :) = locs(eind, :);
	semivari(i+1, :) = semivariogram(data, center(i+1, :), r, h, delta);
	% delete the center and find the second largest one 
	% if it is almost isolated from other data
	notnum = find(isnan(semivari(i+1, :)));
	while (length(semivari(i+1,:)) - nnz(notnum) < 3 )
		rule_new(eind) = eps;
		[extrema, eind] = max(rule_new);
		center(i+1, :) = locs(eind, :);
		semivari(i+1, :) = semivariogram(data, center(i+1, :), r, h, delta);
		notnum = find(isnan(semivari(i+1, :)));
	end

	L(i+1) = LV(2*semivari(i+1, :), h);
	radius(i+1) = 3*L(i+1);
	separation(i+1) = 2/3*radius(i+1);
	i = i+1;
end



label = ones(size(center, 1), 1);



add_center = center(label==1, :);
add_radius = radius(label==1);
add_radius = add_radius(:);

nS = ones(size(grid_all, 1),1);
new_center = [grid_all; add_center];
new_threshold = [threshold; add_radius];

new_nS = [nS; ones(length(add_radius), 1)];



% check basis function matrix
disp('Begin checking basis function matrix')
S = localbasis(locs_PCTM, new_center, new_threshold);



%% delete basis centers that make basis matrix not full column rank
indkeep = size(grid_all, 1);
[~, ind] = licol(S, indkeep);
new_center = new_center(ind,:);
new_threshold = new_threshold(ind);
new_nS = new_nS(ind);


if size(new_center, 1) > 300 | size(new_center, 1) == size(grid_all, 1)
	disp(strcat('The total number of basis functions is: ',...
			 num2str(size(new_center, 1)))); 
	stop_update = 1 ;
else
	stop_update = 0;
end

% calculate residuals with new basis functions 
[res_new] = BPC(locs_PCTM, Y, new_center, new_threshold, sig2eps);



if nargout > 5
	elMSE = rule_old;
end

