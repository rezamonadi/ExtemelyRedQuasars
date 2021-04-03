O=[0,0,0];
A=[0,1,1];
[azimA, elevA, rA] = cart2sph(A(1), A(2), A(3));
azimC = azimA;
elevC= elevA+ deg2rad(5);
rC= rA;
azimD=azimA;
elevD=elevA-deg2rad(5);
rD = rA;
[xC,yC,zC] = sph2cart(azimC, elevC, rC);
[xD,yD,zD] = sph2cart(azimD, elevD, rD);
scatter3([0, A(1), xC, xD], [0, A(2), yC, yD], [0, A(3), zC, zD])