function [mvk,sig2mvk,mvrho,mvsig2]=wrapKrige(index,locs_pred,locs_obs,Yobs,sig2eps,mvsize,nu)

if nargin < 7
	nu = 0.5;
end

center = locs_pred(index,:);
%whitenoise = xi(index);

[mvlocs_obs, mvYobs] = findmv(center, [locs_obs, Yobs], mvsize);

%mvlocs_fine = findmv(center, locs_fine, mvsize);

[rho, sigma2] = MLE_Matern(mvYobs, mvlocs_obs, nu, sig2eps);

[mvk, sig2] = spKrigingMatern(mvYobs, mvlocs_obs, center, nu, rho, sigma2, sig2eps);

if nargout > 1
	sig2mvk=sig2;
end

if nargout > 2
	mvrho = rho;
	mvsig2 = sigma2;
end


end 