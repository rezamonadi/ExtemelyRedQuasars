
clc
% % clear
% %  loading 
X = load('../MedSpec/data.dat');
iw3 = X(:,1);
rew = X(:,2);
kt80 = X(:,3);
X_data = X(:,1:3);
X_normal = normalize(X_data, 'range');
% X_data= X_data;
% %  centering 
MainCenter = median(X_normal);
X_Centered = X_normal - MainCenter;
ERQ = X_normal(((iw3>=4.6) & (rew>=2)),:);
ERQCenter = median(ERQ);
ERQVector = ERQCenter - MainCenter;
[nERQ, c] =size(ERQ);
[nX,c] = size(X);
% building erq vectors and nrmalizing 

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
ind = find(cdf>int32(.80*nERQ));
opening_angle = (bins(ind(1)));
rad2deg((opening_angle))
% in_wedge = ((AllDirections) <= opening_angle);
% in_wedge=in_wedge';
% % % sum(in_wedge)



% % % % % % calculating the density on the grid
ND=500;
x = linspace(min(X_normal(:,1)),max(X_normal(:,1)), ND);
y = linspace(min(X_normal(:,2)),max(X_normal(:,2)), ND);
z = linspace(min(X_normal(:,3)),max(X_normal(:,3)), ND);
[xx,yy,zz] = ndgrid(x,y,z);
xi = [xx(:) yy(:) zz(:)];
Silver = (4/5/nX)^(1/7);
% D = mvksdensity(X_normal, xi, 'Bandwidth', Silver*std(X_normal));
load('D_mesh500_Silver_10.mat')
D_mesh = reshape(D_mesh,size(xx));
% save('D_mesh500_Silver_10.mat', 'D_mesh')
% % ploting the iso-surface on top of data points
% % coeff =[0.85, 0.9,0.95, 1, 1.05, 1.1, 1.15];
% coeff =[ 1];
C =['r', 'b', 'g', 'c', 'm', 'k','y','r'];

isovalue1 =0.5*max(max(max(D_mesh)));
isovalue2 =0.1*max(max(max(D_mesh)));
isovalue3 =0.05*max(max(max(D_mesh)));
isovalue4 =0.01*max(max(max(D_mesh)));
isovalue5 =0.005*max(max(max(D_mesh)));
isovalue6 =0.003*max(max(max(D_mesh)));


surf1 = isosurface(xx,yy,zz, D_mesh, isovalue1);
surf2 = isosurface(xx,yy,zz, D_mesh, isovalue2);
surf3 = isosurface(xx,yy,zz, D_mesh, isovalue3);
surf4 = isosurface(xx,yy,zz, D_mesh, isovalue4);
surf5 = isosurface(xx,yy,zz, D_mesh, isovalue5);
surf6 = isosurface(xx,yy,zz, D_mesh, isovalue6);


% n = cross(ERQVector/norm(ERQVector), [1,0,0]);
% n = cross(ERQVector/norm(ERQVector), [0,1,0]);
n = cross(ERQVector/norm(ERQVector), [0,0,1]);
eps=0.001;
ss = 50;
P0 = MainCenter;
mask = (iw3>=4.6) & (rew>=2);
figure;
scatter3(X_normal(mask,1), X_normal(mask,2), X_normal(mask,3), 10, 'b', 'Marker', 'o')
hold on
s = scatter3(X_normal(~mask,1), X_normal(~mask,2), X_normal(~mask,3), .1, 'k', 'Marker','.')
s.CData = [0.2 0.2 0.2]
s.AlphaData=0.5;
hold on
surf_i = surf1;
Vx = surf_i.vertices(:,1) - P0(1);
Vy = surf_i.vertices(:,2) - P0(2);
Vz = surf_i.vertices(:,3) - P0(3);
mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
surf_new.vertices = surf_i.vertices(mask,:);
surf_new.faces = surf_i.faces(mask,:);
scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), ss, 'Marker', '.');
hold on

surf_i = surf2;
Vx = surf_i.vertices(:,1) - P0(1);
Vy = surf_i.vertices(:,2) - P0(2);
Vz = surf_i.vertices(:,3) - P0(3);
mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
surf_new.vertices = surf_i.vertices(mask,:);
surf_new.faces = surf_i.faces(mask,:);
scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), ss, 'Marker', '.');
hold on

surf_i = surf3;
Vx = surf_i.vertices(:,1) - P0(1);
Vy = surf_i.vertices(:,2) - P0(2);
Vz = surf_i.vertices(:,3) - P0(3);
mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
surf_new.vertices = surf_i.vertices(mask,:);
surf_new.faces = surf_i.faces(mask,:);
scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), ss, 'Marker', '.');
hold on

surf_i = surf4;
Vx = surf_i.vertices(:,1) - P0(1);
Vy = surf_i.vertices(:,2) - P0(2);
Vz = surf_i.vertices(:,3) - P0(3);
mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
surf_new.vertices = surf_i.vertices(mask,:);
surf_new.faces = surf_i.faces(mask,:);
scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), ss, 'Marker', '.');
hold on
surf_i = surf5;
Vx = surf_i.vertices(:,1) - P0(1);
Vy = surf_i.vertices(:,2) - P0(2);
Vz = surf_i.vertices(:,3) - P0(3);
mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
surf_new.vertices = surf_i.vertices(mask,:);
surf_new.faces = surf_i.faces(mask,:);
scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), ss, 'Marker', '.');
hold on
% surf_i = surf6;
% Vx = surf_i.vertices(:,1) - P0(1);
% Vy = surf_i.vertices(:,2) - P0(2);
% Vz = surf_i.vertices(:,3) - P0(3);
% mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
% surf_new.vertices = surf_i.vertices(mask,:);
% surf_new.faces = surf_i.faces(mask,:);
% scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
% hold on

% legend('0.5', '0.1', '0.05', '0.01', '0.005', '0.002')



% view(3); axis tight;
% p =plot3([MainCenter(1), ERQCenter(1)],[MainCenter(2), ERQCenter(2)],... 
% [MainCenter(3), ERQCenter(3)] );
% p.LineWidth=1;
% hold on
% for b=0:8
%     mask = (labels==b);
%     if (b>0)
%         scatter3(X_normal(mask,1), X_normal(mask,2), X_normal(mask,3), 20, c(b,:))
%         hold on
%     else
%         scatter3(X_normal(mask,1), X_normal(mask,2), X_normal(mask,3), 0.5, 'k'  )
%         hold on
%     end
    
% end
% hold  on
lx = sprintf('(i-w3-%.2f)/%.2f', min(iw3), max(iw3)-min(iw3));
ly = sprintf('(rew-%.2f)/%.2f', min(rew), max(rew)-min(rew));
lz = sprintf('(kt80-%.2f)/%.2f', min(kt80), max(kt80)-min(kt80));
set(get(gca, 'XLabel'), 'String', lx);
set(get(gca, 'YLabel'), 'String', ly);
set(get(gca, 'ZLabel'), 'String', lz);
alpha(0.4)
% pp= patch(([4.6, 4.6, 4.6, 4.6]-min(iw3))/(max(iw3)-min(iw3)), [0,0,1,1], [0,1,1,0], [0,0,0,0]);
% hold on
% pp.FaceAlpha=0.1;
% pp.FaceColor='red';
% view(3)
