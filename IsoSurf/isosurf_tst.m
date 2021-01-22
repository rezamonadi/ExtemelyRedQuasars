
clc
% clear
%  loading 
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
in_wedge = ((AllDirections) <= opening_angle);
in_wedge=in_wedge';
% sum(in_wedge)



% % % % % calculating the density on the grid
ND=250;
x = linspace(min(X_normal(:,1)),max(X_normal(:,1)), ND);
y = linspace(min(X_normal(:,2)),max(X_normal(:,2)), ND);
z = linspace(min(X_normal(:,3)),max(X_normal(:,3)), ND);
[xx,yy,zz] = ndgrid(x,y,z);
xi = [xx(:) yy(:) zz(:)];
Silver = (4/5/nX)^(1/7);
% D = mvksdensity(X_normal, xi, 'Bandwidth', 0.8*Silver*std(X_normal));
% D_mesh = reshape(D,size(xx));
% save('D_mesh250_Silver_08.mat', 'D_mesh')
% ploting the iso-surface on top of data points
% coeff =[0.85, 0.9,0.95, 1, 1.05, 1.1, 1.15];
coeff =[ 1];
C =['r', 'b', 'g', 'c', 'm', 'k','y','r'];
for  i=1:1
    isovalue1 =0.5*max(max(max(D_mesh)))*coeff(i);
    isovalue2 =0.1*max(max(max(D_mesh)))*coeff(i);
    isovalue3 =0.05*max(max(max(D_mesh)))*coeff(i);
    isovalue4 =0.01*max(max(max(D_mesh)))*coeff(i);
    surf1 = isosurface(xx,yy,zz, D_mesh, isovalue1);
    surf2 = isosurface(xx,yy,zz, D_mesh, isovalue2);
    surf3 = isosurface(xx,yy,zz, D_mesh, isovalue3);
    surf4 = isosurface(xx,yy,zz, D_mesh, isovalue4);
    % surfe1=surf3;
    % surfe2=surf3;
    % surfe3=surf3;
    % surfe4=surf3;
    % surfe1.vertices=(surf3.vertices - MainCenter)*1.4+ MainCenter;
    % surfe2.vertices=(surf3.vertices - MainCenter)*2.1+ MainCenter;
    % surfe3.vertices=(surf3.vertices - MainCenter)*2.5+ MainCenter;
    % surfe4.vertices=(surf3.vertices - MainCenter)*3+ MainCenter;

    % % initializing labels
    % labels = zeros(1,nX);
    % in_surf1 = inpolyhedron(surf1, X_normal);
    % in_surf2 = inpolyhedron(surf2, X_normal);
    % in_surf3 = inpolyhedron(surf3, X_normal);
    % in_surfe1 = inpolyhedron(surfe1, X_normal);
    % in_surfe2 = inpolyhedron(surfe2, X_normal);
    % in_surfe3 = inpolyhedron(surfe3, X_normal);
    % in_surfe4 = inpolyhedron(surfe4, X_normal);


    % labels((in_surf1==1))=1;
    % labels((in_wedge==1) & (in_surf1==0) & (in_surf2==1))=2;
    % labels((in_wedge==1) & (in_surf2==0) & (in_surf3==1))=3;
    % labels((in_wedge==1) & (in_surf3==0) & (in_surfe1==1))=4;
    % labels((in_wedge==1) & (in_surfe1==0) & (in_surfe2==1))=5;
    % labels((in_wedge==1) & (in_surfe2==0) & (in_surfe3==1))=6;
    % labels((in_wedge==1) & (in_surfe3==0) & (in_surfe4==1))=7;
    % labels((in_wedge==1) & (in_surfe4==0))=8;
    % fid =sprintf('labels-silver-1-open-80-levels-05-01-005-ex-14-21-25-30.mat');
    % save(fid,'labels');
    % c=turbo(8);
    figure;
    % for g=1:6
    %     datax = X_data(labels==g,1);
    %     datay = X_data(labels==g,2);
    %     dataz = X_data(labels==g,3);
    %     % datax = X_normal(labels==g,1);
    %     % datay = X_normal(labels==g,2);
    %     % dataz = X_normal(labels==g,3);
    %     scatter3(datax, datay, dataz,  color(g+1), 'filled')
    %     hold on
    % end
    % % p1 = patch(surf1);
    % % hold on
    % % p1 = patch(surf2);
    % % hold on
    % p1 = patch(surfe3);
    % hold on
    % p1 = patch(surfe3);
    % hold on
    % view(3); axis tight
    % camlight; lighting gouraud
    % alpha(0.3)
    % tit= sprintf('coeff=%.2f', coeff(i));
    % title(tit)


    n = cross(ERQVector/norm(ERQVector), [0,1,0]);
    eps=0.001;
    P0 = MainCenter;

    figure;

    surf_i = surf1;
    Vx = surf_i.vertices(:,1) - P0(1);
    Vy = surf_i.vertices(:,2) - P0(2);
    Vz = surf_i.vertices(:,3) - P0(3);
    mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    surf_new.vertices = surf_i.vertices(mask,:);
    surf_new.faces = surf_i.faces(mask,:);
    scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    hold on

    surf_i = surf2;
    Vx = surf_i.vertices(:,1) - P0(1);
    Vy = surf_i.vertices(:,2) - P0(2);
    Vz = surf_i.vertices(:,3) - P0(3);
    mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    surf_new.vertices = surf_i.vertices(mask,:);
    surf_new.faces = surf_i.faces(mask,:);
    scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    hold on
    
    surf_i = surf3;
    Vx = surf_i.vertices(:,1) - P0(1);
    Vy = surf_i.vertices(:,2) - P0(2);
    Vz = surf_i.vertices(:,3) - P0(3);
    mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    surf_new.vertices = surf_i.vertices(mask,:);
    surf_new.faces = surf_i.faces(mask,:);
    scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    hold on

    surf_i = surf4;
    Vx = surf_i.vertices(:,1) - P0(1);
    Vy = surf_i.vertices(:,2) - P0(2);
    Vz = surf_i.vertices(:,3) - P0(3);
    mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    surf_new.vertices = surf_i.vertices(mask,:);
    surf_new.faces = surf_i.faces(mask,:);
    scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    hold on



    % surf_i = surfe1;
    % Vx = surf_i.vertices(:,1) - MainCenter(1);
    % Vy = surf_i.vertices(:,2) - MainCenter(2);
    % Vz = surf_i.vertices(:,3) - MainCenter(3);
    % mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    % surf_new.vertices = surf_i.vertices(mask,:);
    % surf_new.faces = surf_i.faces(mask,:);
    % scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    % hold on


    % surf_i = surfe2;
    % Vx = surf_i.vertices(:,1) - MainCenter(1);
    % Vy = surf_i.vertices(:,2) - MainCenter(2);
    % Vz = surf_i.vertices(:,3) - MainCenter(3);
    % mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    % surf_new.vertices = surf_i.vertices(mask,:);
    % surf_new.faces = surf_i.faces(mask,:);
    % scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    % hold on


    % surf_i = surfe3;
    % Vx = surf_i.vertices(:,1) - MainCenter(1);
    % Vy = surf_i.vertices(:,2) - MainCenter(2);
    % Vz = surf_i.vertices(:,3) - MainCenter(3);
    % mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    % surf_new.vertices = surf_i.vertices(mask,:);
    % surf_new.faces = surf_i.faces(mask,:);
    % scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    % hold on

    % surf_i = surfe4;
    % Vx = surf_i.vertices(:,1) - MainCenter(1);
    % Vy = surf_i.vertices(:,2) - MainCenter(2);
    % Vz = surf_i.vertices(:,3) - MainCenter(3);
    % mask = abs(Vx*n(1) + Vy*n(2) + Vz*n(3))<eps;
    % surf_new.vertices = surf_i.vertices(mask,:);
    % surf_new.faces = surf_i.faces(mask,:);
    % scatter3(surf_new.vertices(:,1), surf_new.vertices(:,2), surf_new.vertices(:,3), 100, 'Marker', '.');
    % hold on



    % view(3); axis tight;
    p =plot3([MainCenter(1), ERQCenter(1)],[MainCenter(2), ERQCenter(2)], [MainCenter(3), ERQCenter(3)] );
    p.LineWidth=3;
    hold on
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
    % % patch(surf_i)
    lx = sprintf('(i-w3-%.2f)/%.2f', min(iw3), max(iw3)-min(iw3));
    ly = sprintf('(rew-%.2f)/%.2f', min(rew), max(rew)-min(rew));
    lz = sprintf('(kt80-%.2f)/%.2f', min(kt80), max(kt80)-min(kt80));
    set(get(gca, 'XLabel'), 'String', lx);
    set(get(gca, 'YLabel'), 'String', ly);
    set(get(gca, 'ZLabel'), 'String', lz);
    % alpha(0.5)
end
% patch(([4.6, 4.6, 4.6, 4.6]-min(iw3))/(max(iw3)-min(iw3)), [0,0,1,1], [0,1,1,0], [0,0,0,0]);
% xlim([0,1])
% ylim([0,1])
% zlim([0,1])
% sum(labels==1)
% sum(labels==2)
% sum(labels==3)
% sum(labels==4)
% sum(labels==5)
% sum(labels==6)
% sum(labels==7)
% sum(labels==8)