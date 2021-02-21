% Note: This script is used for debug
% Collision detection function (F_CollideCheck) is called under a predefine
% static case. If two ellisoids intersect, the contact point on each
% particle is calculated.
% 
% Ref:Chesnutt & Marshall,Computers & Fluids,38,1782-1794,2009.

clear;
clc;

lplot1 = 1;
lplot2 = 1;

% length of semi axis
a = 2;
b = 1;
c = 1;

% particle centroid: npars*3d
rpar = [3.6,0,0;...
        1,0,0];

% Rotation matrix: 
RMatrx = zeros(2,3,3);
RMatrx(1,:,:) = eye(3); % DO NOT rotate here
RMatrx(2,:,:) = eye(3); % DO NOT rotate here

thickness = 0.0; % added thickness for collision detection
i = 1; % particle id
j = 2; % particle id

% Call function for collision detection and find the contact point
[lcollide,Contact1,Contact2] = F_CollideCheck(i,j,rpar,RMatrx,a,b,c,thickness);

% Plot ellipsoids
nplot = 100;
[xp1,yp1,zp1] = ellipsoid(rpar(i,1),rpar(i,2),rpar(i,3),a,b,c,nplot);
[xp2,yp2,zp2] = ellipsoid(rpar(j,1),rpar(j,2),rpar(j,3),a,b,c,nplot);

figure(1);
if (lplot1)
surf(xp1,yp1,zp1,'FaceAlpha',0.25,'EdgeColor','none');
hold on;
scatter3(Contact1(1),Contact1(2),Contact1(3),10,'r','filled');
axis equal;
end
if (lplot2)
surf(xp2,yp2,zp2,'FaceAlpha',0.25,'EdgeColor','none');
hold on;
scatter3(Contact2(1),Contact2(2),Contact2(3),10,'b','filled');
axis equal;
end

xlabel('X');
ylabel('Y');
zlabel('Z');
