
clc
clear
%  loading 
X = load('data3d.dat');
iw3 = X(:,1);
rew = X(:,2);
kt80 = X(:,3);
X_data = X(:,1:3);
X_normal = normalize(X_data, 'range');
size(X_data);
% X_data= X_data;
% %  centering 
MainCenter = median(X_normal);
X_Centered = X_normal - MainCenter;
ERQ = X_normal(((iw3>=4.6) & (rew>=2) & (kt80>=0.33)),:);
ERQCenter = median(ERQ);
ERQVector = ERQCenter - MainCenter;
[nERQ, c] =size(ERQ);
[nX,c] = size(X);
% % building erq vectors and nrmalizing 
r=.99;
for i=1:nERQ
     ERQNorm = norm(ERQ(i,:));
     ERQDirections(i) = acos(dot(ERQVector, ERQ(i,:))/norm(ERQVector)/ERQNorm);
end
for i=1:nX
    XNorm = norm(X_Centered(i,:));
    AllDirections(i) = acos(dot(ERQVector, X_Centered(i,:))/norm(ERQVector)/XNorm);
end

[counts, bins] = histcounts(ERQDirections, 1000);
cdf = cumsum(counts);
ind = find(cdf>int32(r*nERQ));
opening_angle = (bins(ind(1)));
rad2deg((opening_angle))
in_wedge = ((AllDirections) <= opening_angle);
in_wedge=in_wedge';
% % % % % % % calculating the density on the grid
ND=100;
x = linspace(min(X_normal(:,1)),max(X_normal(:,1)), ND);
y = linspace(min(X_normal(:,2)),max(X_normal(:,2)), ND);
z = linspace(min(X_normal(:,3)),max(X_normal(:,3)), ND);
[xx,yy,zz] = ndgrid(x,y,z);
xi = [xx(:) yy(:) zz(:)];
Silver = (4/5/nX)^(1/7)
std(X_normal);
Silver*std(X_normal);
load('D_mesh100_Silver_10.mat');

    isovalue1 =0.5*max(max(max(D_mesh)));
    isovalue2 =0.05*max(max(max(D_mesh)));
    isovalue3 =0.025*0.65*max(max(max(D_mesh)));

    surf1 = isosurface(xx,yy,zz, D_mesh, isovalue1);
    surf2 = isosurface(xx,yy,zz, D_mesh, isovalue2);
    surf3 = isosurface(xx,yy,zz, D_mesh, isovalue3);
    surfe1=surf3;
    surfe2=surf3;
    surfe3=surf3;
    surfe4=surf3;
    % surfe1.vertices=(surf3.vertices - MainCenter)*1.3+ MainCenter;
    surfe2.vertices=(surf3.vertices - MainCenter)*1.49+ MainCenter;
    surfe3.vertices=(surf3.vertices - MainCenter)*1.95+ MainCenter;
    surfe4.vertices=(surf3.vertices - MainCenter)*2.5+ MainCenter;


t= linspace(0,1,1000);
x_line = ERQVector(1)/norm(ERQVector)*t + MainCenter(1);
y_line = ERQVector(2)/norm(ERQVector)*t + MainCenter(2);
z_line = ERQVector(3)/norm(ERQVector)*t + MainCenter(3);


min_D=100000;
for i=1:1000
    % disp(i)
    for j=1:size(surfe2.vertices(:,1))
        D2 = (surfe2.vertices(j,1)-x_line(i))^2 + (surfe2.vertices(j,2)-y_line(i))^2+...
                (surfe2.vertices(j,3)-z_line(i))^2;
        if(min_D>D2)
            min_D=D2;
            i_min=i;
        end
    end
end
P_intersection = [x_line(i_min), y_line(i_min), z_line(i_min)];
% P_intersection_original = P_intersection.*range(X_data) + min(X_data)
n=ERQVector/norm(ERQVector);
plane = dot(n,P_intersection);
save('Intersection.mat','P_intersection', 'n', 'plane');
