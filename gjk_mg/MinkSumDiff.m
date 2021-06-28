% MINKOWSKI SUM/DIFFERENCE VISUALIZATION
clc; clear all; close all;

l2d = true;

%% BASE GEOMETRIES

% Sphere
% theta = linspace(0,2*pi,20);
% phi = linspace(0,pi,20);
% numpts = length(theta)*length(phi);
% 
% sx1 = zeros(numpts,1);
% sy1 = zeros(numpts,1);
% sz1 = zeros(numpts,1);
% 
% for t=1:length(theta)
%     for p=1:length(phi)
%         th = theta(t);
%         ph = phi(p);
%         
%         sx1((t-1)*length(theta) + p) = cos(th)*sin(ph);
%         sy1((t-1)*length(theta) + p) = sin(th)*sin(ph);
%         sz1((t-1)*length(theta) + p) = cos(ph);
%     end
% end
% 
% sx1 = 2.5 + sx1;
% sy1 = 2.5 + sy1;
% sz1 = 18 + sz1;

% Circle
% theta = linspace(0,2*pi,150);
% sx1 = 2 + cos(theta); sy1 = sin(theta); sz1 = zeros(length(theta),1)';

% Triangle
%Shape1 = [0 1 2 3 4 5; 0 2.5 5 2.5 0; 0 2.5 5 2.5 0];
%Shape1 = [linspace(0,5,50), linspace(5,0,25); linspace(0,5,25), linspace(5,0,25), zeros(25,1)'; linspace(0,5,25), linspace(5,0,25), zeros(25,1)'];

% Pyramid
% sx1 = [-1,1,0,0]; sy1 = [-1,-1,1,0]; sz1 = [0,0,0,1];
% dx = 0; dy = 0; dz = 0;
% sx1 = sx1+dx; sy1 = sy1+dy; sz1 = sz1+dz; 

% Square
% sx1 = []; sy1 = []; 



%% Shape 1

theta = linspace(0,2*pi,150);
sx1 = 0.1+cos(theta); sy1 = 0.2+ sin(theta); sz1 = zeros(length(theta),1)';

numpts1 = length(sx1);

%% Shape 2

%Shape2 = [linspace(0,5,50), linspace(5,0,25); linspace(0,5,25), linspace(5,0,25), zeros(25,1)'; linspace(0,5,25), linspace(5,0,25), zeros(25,1)'];

Shape2 = [-4,-4,0,0; -4,0,0,-4; 0,0,0,0];

sx2 = Shape2(1,:); sy2 = Shape2(2,:); sz2 = Shape2(3,:);
numpts2 = length(sx2);

%% Shape Visualization
figure;
hold on;
grid on;
box on;
axis equal;
set(gcf,'position',[0,400,400,400]);

if(l2d==false)
    k1 = convhull(sx1,sy1,sz1);
    trisurf(k1,sx1,sy1,sz1);
    k2 = convhull(sx2,sy2,sz2);
    trisurf(k2,sx2,sy2,sz2);
else
    patch(sx1,sy1,sz1,'r');
    patch(sx2,sy2,sz2,'b');
end


view(45,45);


%% Minkowski Sum Visualization

figure;
hold on;
grid on;
box on;
axis equal;
set(gcf,'position',[400,400,400,400]);
view(45,30);

if(l2d==false)
    k1 = convhull(sx1,sy1,sz1);
    trisurf(k1,sx1,sy1,sz1);
    k2 = convhull(sx2,sy2,sz2);
    trisurf(k2,sx2,sy2,sz2);
else
    patch(sx1,sy1,sz1,'r');
    patch(sx2,sy2,sz2,'b');
end

numptsSum = numpts1*numpts2;
msumx = zeros(numptsSum,1); msumy = zeros(numptsSum,1); msumz = zeros(numptsSum,1);

% %Shape1 - Shape2
% for i = 1:numpts1
%     for j = 1:numpts2
%        sdx_ij = sx1(i) - sx2(j);
%        sdy_ij = sy1(i) - sy2(j);
%        sdz_ij = sz1(i) - sz2(j);
%        
%        sdx((i-1)*numpts1+j) = sdx_ij;
%        sdy((i-1)*numpts1+j) = sdy_ij;
%        sdz((i-1)*numpts1+j) = sdz_ij;
%     end
% end


%Shape2 - Shape1
for i = 1:numpts2
    for j = 1:numpts1
       sdx_ij = sx2(i) + sx1(j);
       sdy_ij = sy2(i) + sy1(j);
       sdz_ij = sz2(i) + sz1(j);
       
       msumx((i-1)*numpts1+j) = sdx_ij;
       msumy((i-1)*numpts1+j) = sdy_ij;
       msumz((i-1)*numpts1+j) = sdz_ij;
    end
end


if(l2d==false)
    ksum = convhull(msumx,msumy,msumz);
    trimesh(ksum,msumx,msumy,msumz);
else
    scatter3(msumx,msumy,msumz,'mp');
end

title('\it{Minkowski Sum}','fontsize',16,'interpreter','latex');




%% Minkowski Difference Visualization

figure;
hold on;
grid on;
box on;
axis equal;
set(gcf,'position',[0,0,400,400]);

if(l2d==false)
    k1 = convhull(sx1,sy1,sz1);
    trisurf(k1,sx1,sy1,sz1);
    k2 = convhull(sx2,sy2,sz2);
    trisurf(k2,sx2,sy2,sz2);
else
    patch(sx1,sy1,sz1,'r');
    patch(sx2,sy2,sz2,'b');
end


numptsDiff = numpts1*numpts2;
mdiffx = zeros(numptsDiff,1); mdiffy = zeros(numptsDiff,1); mdiffz = zeros(numptsDiff,1);

% %Shape1 - Shape2
% for i = 1:numpts1
%     for j = 1:numpts2
%        sdx_ij = sx1(i) - sx2(j);
%        sdy_ij = sy1(i) - sy2(j);
%        sdz_ij = sz1(i) - sz2(j);
%        
%        sdx((i-1)*numpts1+j) = sdx_ij;
%        sdy((i-1)*numpts1+j) = sdy_ij;
%        sdz((i-1)*numpts1+j) = sdz_ij;
%     end
% end


%Shape2 - Shape1
for i = 1:numpts2
    for j = 1:numpts1
       sdx_ij = sx2(i) - sx1(j);
       sdy_ij = sy2(i) - sy1(j);
       sdz_ij = sz2(i) - sz1(j);
       
       mdiffx((i-1)*numpts1+j) = sdx_ij;
       mdiffy((i-1)*numpts1+j) = sdy_ij;
       mdiffz((i-1)*numpts1+j) = sdz_ij;
    end
end

%scatter3(sdx,sdy,sdz,'filled');

% patch(sdx, sdy, sdz,'c');
if(l2d==false)
    kdiff = convhull(mdiffx,mdiffy,mdiffz);
    trimesh(kdiff,mdiffx,mdiffy,mdiffz);
else
    scatter3(mdiffx,mdiffy,mdiffz,50,'m','filled');
end

view(45,30);
set(gcf,'position',[400,0,400,400]);
title('\it{Minkowski Difference}','fontsize',16,'interpreter','latex');


