X = load('../MedSpec/data.dat');
iw3 = X(:,1);
rew = X(:,2);
[nX,c] = size(X);

X_normal = normalize(X(:,1:2), 'range');
ND=250;
x = linspace(min(X_normal(:,1)),max(X_normal(:,1)), ND);
y = linspace(min(X_normal(:,2)),max(X_normal(:,2)), ND);
[xx,yy] = ndgrid(x,y);
xi = [xx(:) yy(:)];
Silver = (1/nX)^(1/6)
std(X_normal)
Silver*std(X_normal)
% D = mvksdensity(X_normal, xi, 'Bandwidth', Silver*std(X_normal));
% % load('2d_D_mesh250_Silver_08.mat')
% D_mesh = reshape(D,size(xx));
% D_mesh=D_mesh/(max(max(D_mesh)));
% save('2d_D_mesh250_Silver_08.mat', 'D_mesh')
% mask = (iw3>=4.6) & (rew>=2);
% scatter(X_normal(mask,1), X_normal(mask,2), 10, 'b', 'Marker', 'o')
% hold on
% scatter(X_normal(~mask,1), X_normal(~mask,2), 1, 'k', 'Marker', '.')
% hold on

% [M,c] = contour(xx,yy,D_mesh, [0.5, 0.1,0.05, 0.01,0.005, 0.0015]);
% clabel(M,c, 'manual',  'FontSize',10, 'Color', 'r')
% c.LineWidth=2;
% colormap('turbo')
% lx = sprintf('(i-w3-%.2f)/%.2f', min(iw3), max(iw3)-min(iw3));
% ly = sprintf('(rew-%.2f)/%.2f', min(rew), max(rew)-min(rew));
% set(get(gca, 'XLabel'), 'String', lx);
% set(get(gca, 'YLabel'), 'String', ly);
% alpha(0.7)
