	
% Author: Pulong Ma <mpulong@gmail.com>
% Date: March 05, 2017 
% Last Modified by: Pulong Ma
% Last Modified Date: March 05, 2017
% Last Modified time: 12:16:24


% Purpose: find linearly independent columns after first few columns




function [Xsub, ind] = licol(S, indkeep, tol);
%%% 
% S: an n-by-r matrix, with r < n.
% indkeep: an integer indicating keeping the first few columns.

if nargin < 3
	tol=1e-4;
end

[n r] = size(S);
ind = 1:indkeep;
count = 1;
for i=(indkeep+1):r
	temp = S(:,[ind,i]).'*S(:,[ind,i]);
	if rank(full(temp), tol)==(indkeep+count)
		ind(indkeep+count) = i;
		count = count + 1;
	end
end
ind = ind(:);
Xsub = S(:,ind);



