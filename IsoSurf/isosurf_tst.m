
clc
clear
%  loading 
X = load('../MedSpec/data.dat');
iw3 = X(:,1);
rew = X(:,2);
kt80 = X(:,3);
X_data = X(:,1:3);
X_normal = normalize(X_data);
% %  centering 
MainCenter = median(X_normal);
X_Centered = X_normal - MainCenter;
ERQ = X_normal(((iw3>=4.6) & (rew>=2)),:);
ERQCenter = median(ERQ);
ERQVector = ERQCenter -MainCenter;
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
in_wedge = ((AllDirections) <= opening_angle);
in_wedge=in_wedge';
sum(in_wedge)



% % % % calculating the density on the grid
ND=25;
x = linspace(min(iw3),max(iw3), ND);
y = linspace(min(rew),max(rew), ND);
z = linspace(min(kt80),max(kt80), ND);
[xx,yy,zz] = ndgrid(x,y,z);
xi = [xx(:) yy(:) zz(:)];
D = mvksdensity(X_data, xi, 'Bandwidth',[0.6*std(iw3),0.6*std(rew),0.6*std(kt80)]);
D_mesh = reshape(D,size(xx));

% ploting the iso-surface on top of data points
% figure
MainCenter = median(X_data);
coeef =[0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3]
for  i=1:7
    isovalue1 =1/3*max(D)*coeef(i);
    isovalue2 =1/9*max(D)*coeef(i);
    surf1 = isosurface(xx,yy,zz, D_mesh, isovalue1);
    surf2 = isosurface(xx,yy,zz, D_mesh, isovalue2);
    surfe1=surf2;
    surfe2=surf2;
    surfe1.vertices=(surfe1.vertices - MainCenter)*1.5+ MainCenter;
    surfe2.vertices=(surfe2.vertices - MainCenter)*2+ MainCenter;

    % initializing labels
    labels = zeros(1,nX);
    in_surf1 = inpolyhedron(surf1, X_data);
    in_surf2 = inpolyhedron(surf2, X_data);
    in_surfe1 = inpolyhedron(surfe1, X_data);
    in_surfe2 = inpolyhedron(surfe2, X_data);
    labels((in_surf1==1))=1;
    labels((in_wedge==1) & (in_surf1==0) & (in_surf2==1))=2;
    labels((in_wedge==1) & (in_surf2==0) & (in_surfe1==1))=3;
    labels((in_wedge==1) & (in_surfe1==0) & (in_surfe2==1))=4;
    labels((in_wedge==1) & (in_surfe2==0))=5;
    fid =sprintf('labels-coeff%.1f-inv3-inv9-exp1-15-exp2-2-open-80.mat', coeef(i))
    save(fid,'labels')
end
% color = ['c','y','c','m','b','r'];
% for g=1:5
%     datax = X_data(labels==g,1);
%     datay = X_data(labels==g,2);
%     dataz = X_data(labels==g,3);
%     % datax = X_normal(labels==g,1);
%     % datay = X_normal(labels==g,2);
%     % dataz = X_normal(labels==g,3);
%     scatter3(datax, datay, dataz,  color(g+1), 'filled')
%     hold on
% end
% p1 = patch(surf1);
% hold on
% p1 = patch(surf2);
% hold on
% p1 = patch(surfe1);
% hold on
% p1 = patch(surfe2);
% hold on
% % % % % % p1 = patch(surfe2);
% % % % % % hold on
% % % % % % alpha(0.1)
% % % scatter3(X_normal(:,1), X_normal(:,2),  X_normal(:,3), 5, 'Marker', '.')
% % % scatter3(MainCenter(1), MainCenter(2), MainCenter(3), 1000, 'Marker', '.')
% % % % daspect([1,1,1])
% view(3); axis tight
% camlight; lighting gouraud
% alpha(0.3)
