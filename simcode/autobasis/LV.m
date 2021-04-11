% Author: Pulong Ma <mpulong@gmail.com>
% Date: January 20th, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: January 20th, 2015
% Last Modified time: 23:00:00


% Purpose: 



function [range, sig2] = LV(vari, h)
% input: variogram and lag distance


range = zeros(size(vari, 1), 1);
sig2 = zeros(size(vari, 1), 1);
%lb = [0.0001; 0.165];
lb = [0.0001; 0.002];
for i = 1:size(vari, 1)
	xdata = h';
	ydata = vari(i, :)'/2;
	k = find(isnan(ydata));
	xdata(k) = [];
	ydata(k) = [];
	x0 = [max(ydata)/2; 3];
	ub = [max(ydata)*2; 12/3];
	myfun = @(x, xdata) x(1)*(1 - exp(-xdata./x(2)));
	options = optimset('MaxIter', 1000, 'TolFun', 1e-4, 'Display', 'off');
	[x, resnorm] = lsqcurvefit(myfun, x0, xdata, ydata, lb, ub, options);
	sig2(i) = x(1);	
	range(i) = x(2);
end

if nargout > 1
	sigma2 = sig2;
end
end