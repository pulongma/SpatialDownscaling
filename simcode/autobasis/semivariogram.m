% Author: Pulong Ma <mpulong@gmail.com>
% Date: March 2nd, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: April 25th, 2015
% Last Modified time: 22:57:11


% Purpose: calculate traditional empirical semivariogram or
% robust semivariogram 




function [semivari, n_pairs] = semivariogram(data, center, r, h, delta, type)

%% Input Arguments:
% (1) data: a cell array with first element being lattitude, second 
% element being longitude, third element being residuals of the data
% (2) center: centers of region where the semivariogram is calculated.
% (3) r: range/radius the the region
% (4) h: bin size to calculate semivariogram.
% (5) type: 
% type=0: traditional empirical semivariogram estimator
% type=1: robust estimator for semivariogram in Cressie and Hawkins(1980)

if nargin < 3
	r = 1000;
	h = 50:50:2000;
	delta = 25;
	type = 0;
elseif nargin < 6
	type = 1;
end	

lat = data{1};
lon = data{2};
locs = [lat, lon];
residuals = data{3};

N = zeros(size(center, 1), length(h));
semivari = N;

switch(type)
case 0
for i = 1:size(center, 1)
	dtemp = great_circle_distance([lat, lon], center(i, :));
	dtemp(locs(:,1)==center(i,1) & locs(:,2)==center(i,2)) = 0;
	indtemp = find(dtemp <= r);
	res = residuals(indtemp);
	locs_temp = locs(indtemp, :);
	dist = zeros(length(res), length(res));
	for j = 1:length(res)
		d = great_circle_distance(locs_temp, locs_temp(j, :));
		tag = find(locs_temp(:,1)==locs_temp(j,1) & locs_temp(:,2)==locs_temp(j,2));
		d(tag) = 0;
		dist(j:end, j) = d(j:end);  % lower triangular distance matrix
	end

	for k = 1:length(h)
		[indi, indj] = find(dist > h(k)-delta & dist <= h(k)+delta);
		N(i, k) = length(indi);
		semivari(i, k) = mean((res(indi) - res(indj)).^2)/2;
	end
	clear dtemp indtemp res locs_temp dist d tag;
end

case 1
for i = 1:size(center, 1)
	dtemp = great_circle_distance([lat, lon], center(i, :));
	dtemp(locs(:,1)==center(i,1) & locs(:,2)==center(i,2)) = 0;
	indtemp = find(dtemp <= r);
	res = residuals(indtemp);
	locs_temp = locs(indtemp, :);
	dist = zeros(length(res), length(res));
	for j = 1:length(res)
		d = great_circle_distance(locs_temp, locs_temp(j, :));
		tag = find(locs_temp(:,1)==locs_temp(j,1) & locs_temp(:,2)==locs_temp(j,2));
		d(tag) = 0;
		dist(j:end, j) = d(j:end);  % lower triangular distance matrix
	end

	for k = 1:length(h)
		[indi, indj] = find(dist > h(k)-delta & dist <= h(k)+delta);
		N(i, k) = length(indi);
		semivari(i, k) = 0.5*(mean(abs(res(indi) - res(indj)).^0.5)).^4 ...
						 ./ (0.457 + 0.494/N(i, k));
	end
	clear dtemp indtemp res locs_temp dist d tag;
end

end

if nargout > 1
	n_pairs = N;
end

