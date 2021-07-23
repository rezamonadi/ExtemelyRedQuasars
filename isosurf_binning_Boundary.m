
clc
clear
%  loading 
X = load('data3d.dat');
iw3 = X(:,1);
rew = X(:,2);
kt80 = X(:,3);
X_data = X(:,1:3);
X_normal = normalize(X_data, 'range');
size(X_data)
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
% sum(in_wedge)



% % % % % % % calculating the density on the grid
ND=100;
x = linspace(min(X_normal(:,1)),max(X_normal(:,1)), ND);
y = linspace(min(X_normal(:,2)),max(X_normal(:,2)), ND);
z = linspace(min(X_normal(:,3)),max(X_normal(:,3)), ND);
[xx,yy,zz] = ndgrid(x,y,z);
xi = [xx(:) yy(:) zz(:)];
Silver = (4/5/nX)^(1/7)
std(X_normal)
Silver*std(X_normal)
% % D = mvksdensity(X_normal, xi, 'Bandwidth', Silver*std(X_normal));
load('D_mesh100_Silver_10.mat');
% % D_mesh = reshape(D, size(xx));
% % save('D_mesh100_Silver_10.mat', 'D_mesh')

% ploting the iso-surface on top of data points

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
surfe4.vertices=(surf3.vertices - MainCenter)*2.5+ MainCenter;
for ii=-10:10
    disp(ii)
    for jj=-10:10
        surfe2.vertices=(surf3.vertices - MainCenter)*(1.5+ii*0.01)+ MainCenter;
        surfe3.vertices=(surf3.vertices - MainCenter)*(2.0+jj*0.01)+ MainCenter;

        % initializing labels
        labels = zeros(1,nX);
        in_surf1 = inpolyhedron(surf1, X_normal);
        in_surf2 = inpolyhedron(surf2, X_normal);
        in_surf3 = inpolyhedron(surf3, X_normal);
        % in_surfe1 = inpolyhedron(surfe1, X_normal);
        in_surfe2 = inpolyhedron(surfe2, X_normal);
        in_surfe3 = inpolyhedron(surfe3, X_normal);
        in_surfe4 = inpolyhedron(surfe4, X_normal);

        labels((in_surf1==1))=1;
        labels((in_wedge==1) & (in_surf1==0) & (in_surf2==1))=2;
        labels((in_wedge==1) & (in_surf2==0) & (in_surf3==1))=3;
        labels((in_wedge==1) & (in_surf3==0) & (in_surfe2==1))=4;
        labels((in_wedge==1) & (in_surfe2==0) & (in_surfe3==1))=5;
        % labels((in_wedge==1) & (in_surfe2==0) & (in_surfe3==1))=6;
        labels((in_wedge==1) & (in_surfe3==0) & (in_surfe4==1))=6;
        labels((in_wedge==1) & (in_surfe4==0))=7;
        save(sprintf('labels-%d-%d.mat',ii,jj),'labels')

    end
end
