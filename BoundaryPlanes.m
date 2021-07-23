% clc
clear
clf
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
MainCenter_original = median(X_data);
ERQ_original = X_data(((iw3>=4.6) & (rew>=2) & (kt80>=0.33)),:);
ERQCenter_original = median(ERQ_original);
ERQVector_original = ERQCenter_original -MainCenter_original;
MainCenter = median(X_normal);
X_Centered = X_normal - MainCenter;
ERQ = X_normal(((iw3>=4.6) & (rew>=2)),:);
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
rad2deg((opening_angle));
in_wedge = ((AllDirections) <= opening_angle);
in_wedge=in_wedge';
% -------------------------

% the coordinate system is centered at the center of T1LM
% Getting azimuth, elevation, and r of ERQVector in normalized space
[azimERQVector, ...
    elevERQVector, ...
    rERQVector] = cart2sph(ERQVector(1),...
    ERQVector(2), ERQVector(3));

% Coordinate of the point with higher eleveation than ERQCenter
% in Cartisian system and in the normalized space
[P_normalized_cart_elev_up(1), P_normalized_cart_elev_up(2),...
    P_normalized_cart_elev_up(3)]...
    = sph2cart(azimERQVector,...
    elevERQVector+opening_angle,...
    1.6*rERQVector);

% Getting the vector connecting ERQCenter to the 
% higher elevation point in the normalized space

P_normalized_cart_elev_up_0 = P_normalized_cart_elev_up + MainCenter;
vec_up = P_normalized_cart_elev_up_0 - ERQCenter;
% Getting left vector from vec_up
vec_left = cross(vec_up,ERQVector);
vec_left = vec_left/norm(vec_left);
vec_left = vec_left*norm(vec_up);

% Coordinate of the point with lower eleveation than ERQCenter
% in Cartisian system and in the normalized space
[P_normalized_cart_elev_down(1),...
    P_normalized_cart_elev_down(2),...
    P_normalized_cart_elev_down(3)]...
    = sph2cart(azimERQVector,...
    elevERQVector-opening_angle,...
    1.6*rERQVector);
% vecotr connecting ERQCenter to the point with lower elevation 
% in the cartisien normalized space

P_normalized_cart_elev_down_0 = P_normalized_cart_elev_down+ MainCenter;
vec_down = P_normalized_cart_elev_down_0 - ERQCenter;
% Getting right vector from vec_down
vec_right = cross(vec_down,ERQVector);
vec_right = vec_right/norm(vec_right);
vec_right = vec_right*norm(vec_down);

% Getting a point left of ERQCeneter on the cross-section of cone
P_normalized_cart_left = ERQCenter + vec_left;

% Getting a point right of ERQCeneter on the cross-section of cone
P_normalized_cart_right = ERQCenter + vec_right;

maskERQ = (iw3>=4.6) & (rew>=2);
h=scatter3(X_normal(maskERQ,1), X_normal(maskERQ,2), X_normal(maskERQ,3),20);
hold on 
scatter3(ERQCenter(1), ERQCenter(2),...
    ERQCenter(3), 100)
hold on 
scatter3(P_normalized_cart_elev_down_0(1), ...
    P_normalized_cart_elev_down_0(2),...
    P_normalized_cart_elev_down_0(3))
hold on 
scatter3(P_normalized_cart_elev_up_0(1), ...
    P_normalized_cart_elev_up_0(2),...
    P_normalized_cart_elev_up_0(3))
hold on 

scatter3(P_normalized_cart_left(1), ...
    P_normalized_cart_left(2),...
    P_normalized_cart_left(3))
hold on 
scatter3(P_normalized_cart_right(1), ...
    P_normalized_cart_right(2),...
    P_normalized_cart_right(3))
hold on 

plot3([MainCenter(1), ERQCenter(1)], [MainCenter(2), ERQCenter(2)], [MainCenter(3), ERQCenter(3)])
hold on 



XX = [MainCenter(1), P_normalized_cart_elev_down_0(1), P_normalized_cart_left(1)];
YY = [MainCenter(2), P_normalized_cart_elev_down_0(2), P_normalized_cart_left(2)];
ZZ = [MainCenter(3), P_normalized_cart_elev_down_0(3), P_normalized_cart_left(3)];
pp=patch(XX,YY,ZZ,[.5,0.5,.9]);
pp.FaceAlpha=0.2;

hold on 


XX = [MainCenter(1), P_normalized_cart_elev_down_0(1), P_normalized_cart_right(1)];
YY = [MainCenter(2), P_normalized_cart_elev_down_0(2), P_normalized_cart_right(2)];
ZZ = [MainCenter(3), P_normalized_cart_elev_down_0(3), P_normalized_cart_right(3)];
pp=patch(XX,YY,ZZ,[.5,0.5,.9]);
pp.FaceAlpha=0.2;

hold on 

XX = [MainCenter(1), P_normalized_cart_elev_up_0(1), P_normalized_cart_right(1)];
YY = [MainCenter(2), P_normalized_cart_elev_up_0(2), P_normalized_cart_right(2)];
ZZ = [MainCenter(3), P_normalized_cart_elev_up_0(3), P_normalized_cart_right(3)];
pp=patch(XX,YY,ZZ,[.5,0.5,.9]);
pp.FaceAlpha=0.2;

hold on 
XX = [MainCenter(1), P_normalized_cart_elev_up_0(1), P_normalized_cart_left(1)];
YY = [MainCenter(2), P_normalized_cart_elev_up_0(2), P_normalized_cart_left(2)];
ZZ = [MainCenter(3), P_normalized_cart_elev_up_0(3), P_normalized_cart_left(3)];
pp=patch(XX,YY,ZZ,[.5,0.5,.9]);
pp.FaceAlpha=0.2;

hold on 

%  intersection of the ERQVector and the outer
%  iso-surface of Bin-4 in the normalized space
P = [ 0.5083    0.7271    0.7368];
scatter3(P(1), P(2), P(3), 100)

hold on
v1 = cross(ERQVector, [1,0,0]);
v2 = cross(ERQVector, [0,1,1]);
v3 = cross(ERQVector, [-1,0,0]);
v4 = cross(ERQVector, [0,-1,-1]);
G1 = P + v1*2;
G2 = P+ v2*3;
G3 = P+ v3*2;
G4 = P+ v4*3;

pp= scatter3([G1(1), G2(1), G3(1), G4(1)], [G1(2), G2(2), G3(2), G4(2)] ,[G1(3), G2(3)...
    G3(3), G4(3)]);
pp.Marker= 'x';

XX = [G1(1), G2(1), G3(1), G4(1)];
YY = [G1(2), G2(2), G3(2), G4(2)];
ZZ = [G1(3), G2(3), G3(3), G4(3)];

% pp=patch(XX,YY,ZZ,[1,0,1]);
% pp.FaceAlpha=0.1;



% axes equal




P_orginal = P.*(max(X_data) -min(X_data)) + min(X_data);
n= ERQVector_original/norm(ERQVector_original);

P_up_org = P_normalized_cart_elev_up_0.*(max(X_data) -min(X_data)) + min(X_data);
vec_up_org = vec_up.*(max(X_data) -min(X_data)) + min(X_data);
n_up = vec_up_org/norm(vec_up_org);

P_down_org = P_normalized_cart_elev_down_0.*(max(X_data) -min(X_data)) + min(X_data);
vec_down_org = vec_down.*(max(X_data) -min(X_data)) + min(X_data);
n_down = vec_down_org/norm(vec_down_org);

P_right_org = P_normalized_cart_right.*(max(X_data) -min(X_data)) + min(X_data);
vec_right_org = vec_right.*(max(X_data) -min(X_data)) + min(X_data);
n_right = vec_right_org/norm(vec_right_org);

P_left_org = P_normalized_cart_left.*(max(X_data) -min(X_data)) + min(X_data);
vec_left_org = vec_left.*(max(X_data) -min(X_data)) + min(X_data);
n_left = vec_left_org/norm(vec_left_org);



