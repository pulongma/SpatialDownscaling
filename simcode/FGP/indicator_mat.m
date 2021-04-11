function mat = indicator_mat(site_obs, site)

n = size(site_obs, 1);
N = size(site, 1);
mat = spalloc(n, N, N);

[LIA, LOCB] = ismember(site_obs, site, 'rows', 'legacy');

for k=1:n
	mat(k, LOCB(k)) = 1;
end

% isequal(siteobs, mat*site) % sanity check
