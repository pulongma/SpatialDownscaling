%%%  variogram-estimation of the measurement-error
% (and fine-scale variation)

function [sigma2_epsilon, sigma2_delta, variogram_plot_data]=...
    variogram_estimation(z_tilde,lon,lat,bins,V_eps)
% z_tilde is the (detrended) data
% bins are the endpoints for the bins used for variogram estimation

% use default (identity matrix) for measurement error matrix if none given
if nargin<5,
    v_eps=ones(length(z),1);
else
    v_eps=diag(V_eps);
end

% declare new function that has the right output for pdist()
function distkm=distfun_km(x,y)
    distkm=(distance_spherical(x,y))';
end

% calculate distances between measurement locations
d=round(pdist([lon lat],@distfun_km)); 
D=squareform(d);

% create robust estimator of semivariogram
rob_meanb=[]; n_pairsb=[];
for i=1:(length(bins)-1),
    [I,J] = find(D>bins(i) & D<bins(i+1));
    rob_meanb(i)=.5*mean(sqrt( abs(z_tilde(I)./sqrt(v_eps(I))-...
        z_tilde(J)./sqrt(v_eps(J))) ))^4 / (.457+ .494/length(I));
    n_pairsb(i)=length(I);
end    

% fit straight line to semivariogram using WLS
centers=mean([bins(1:(end-1))' bins(2:end)'],2)';
X_variogram=[ones(length(centers),1) centers'];
beta_variogram=(X_variogram'*diag(n_pairsb)*X_variogram)\...
    (X_variogram'*diag(n_pairsb)*rob_meanb');    

% the measurement-error variance is equal to the intercept of the WLS line
sigma2_epsilon=beta_variogram(1);

% now the variance of the fine-scale variation
[I,J] = find(D>bins(1) & D<bins(2));
sigma2_delta=.5*mean( (z_tilde(I)-z_tilde(J)).^2 - ...
    sigma2_epsilon*(v_eps(I)+v_eps(J)) );

% export plot-data for plotting later
variogram_plot_data={centers,rob_meanb,beta_variogram,n_pairsb,bins};

end