function dist=great_circle_distance(coor1, coor2)
%[xlat_sub(ni) ylon_sub(ni)],[xlat_sub(nj) ylon_sub(nj)]);

%coor1=[xlat_sub(ni) ylon_sub(ni)];
%coor2=[xlat_sub(nj) ylon_sub(nj)];

xlat1=coor1(:,1);
ylon1=coor1(:,2);

xlat2=coor2(:,1);
ylon2=coor2(:,2);

dist = sqrt((xlat1-xlat2).^2+(ylon1-ylon2).^2);