% Author: Pulong Ma <mpulong@gmail.com>
% Date: Febuarary 13th, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: April 5th, 2015
% Last Modified time: 16:44:50


% Purpose: 



function S = localbasis(locs, basis_center, threshold)
n = size(locs, 1);	
r = size(basis_center, 1);
S = sparse(n, r);

for j = 1:r
	d = great_circle_distance(locs, basis_center(j, :));
	[indi, indj] = find(locs(:, 1) == basis_center(j, 1) & ...
					locs(:, 2) == basis_center(j, 2));
	d(indi) = 0;
	S(:, j) = ((1-d./threshold(j)).^2).^2.*(d<=threshold(j));
end
