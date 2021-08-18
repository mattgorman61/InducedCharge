clc; clear all; close all;

dx = [0 0 0 0 0 0]; dy = [0 0 2.1 -0.4 0 2.05 ]; dz = [0 0 0 0 0 0.2];

load('ParData_B0p2_beta6_rotateX_2884.mat');
x1 = p(:,1) + dx(1); 
y1 = p(:,2) + dy(1); 
z1 = p(:,3) + dz(1);

load('ParData_B0p2_beta6_rotateZ_NegPiDiv12_21924.mat');
x2 = p(:,1) + dx(2); 
y2 = p(:,2) + dy(2); 
z2 = p(:,3) + dz(2);

load('ParData_B0p2_beta6_rotateZ_piDiv12_21924.mat');
x3 = p(:,1) + dx(3); 
y3 = p(:,2) + dy(3); 
z3 = p(:,3) + dz(3);

load('ParData_B0p2_beta6_rotateZ_piDiv12_21924.mat');
x4 = p(:,1) + dx(4); 
y4 = p(:,2) + dy(4); 
z4 = p(:,3) + dz(4);

load('ParData_A0p25_alpha15_18854.mat');
x5 = p(:,1) + dx(5); 
y5 = p(:,2) + dy(5); 
z5 = p(:,3) + dz(5);

load('ParData_A0p25_alpha15_18854.mat');
x6 = p(:,1) + dx(6); 
y6 = p(:,2) + dy(6); 
z6 = p(:,3) + dz(6);

% 
% MV - Theta Symm
figure;
axis equal; hold on; grid on; box on;
scatter3(x3,y3,z3,50,'b','filled');
scatter3(x2,y2,z2,50,'r','filled');
view(90,90);

% TT - Theta Symm
figure;
axis equal; hold on; grid on; box on;
scatter3(x3,y3,z3,50,'b','filled');
scatter3(x4,y4,z4,50,[.3 .3 .3],'filled');
view(90,90);

% % MV - Phi Symm
% figure;
% axis equal; hold on; grid on; box on;
% scatter3(x5,y5,z5,10,'b','filled');
% scatter3(x6,y6,z6,10,[.3 .3 .3],'filled');
% view(90,0);





%{
figure;
scatter3(x_save_mg,y_save_mg,z_save_mg, 10, sigma_save_mg);
box on; grid on; axis equal;
view(10,10);

figure;
scatter3(x_save_mg(1:length(x_save_mg)/2),y_save_mg(1:length(y_save_mg)/2),z_save_mg(1:length(z_save_mg)/2),...
    10, sigma_save_mg(1:length(sigma_save_mg)/2));
box on; grid on; axis equal;
view(10,10);

figure;
scatter3(x_save_mg(length(x_save_mg)/2+1:length(x_save_mg)),y_save_mg(length(y_save_mg)/2+1:length(y_save_mg)),z_save_mg(length(z_save_mg)/2+1:length(z_save_mg)),...
    10, sigma_save_mg(length(sigma_save_mg)/2+1:length(sigma_save_mg)));
box on; grid on; axis equal;
view(10,10);

%}