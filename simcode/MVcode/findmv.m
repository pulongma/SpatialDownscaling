function [locs,obs]=findmv(center,data,mvsize)


center_lat = center(1,1);
center_lon = center(1,2);

lat = data(:,1);
lon = data(:,2);


delta = mvsize;

flag = (lat >= center_lat-delta) & (lat <= center_lat+delta) & ...
       (lon >= center_lon-delta) & (lon <= center_lon+delta);
%flag = sqrt((lat-center_lat).^2+(lon-center_lon).^2) <= delta;

locs = [lat(flag), lon(flag)];

if nargout > 1
	obs = data(flag,3);
end


end
