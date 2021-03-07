close all;
clear all;
clc;

% Definition
sigma_f = 0; % Surface free charge density
kappa_air = 1; % Dielectric constance of air
kappa_p = 2.5; % Dielectric constance of target
epsilon = 1; %Unit: F/m

% Geometry: sphere + point charge
Radius = 5; % Sphere radius: 100 micron
N = 5005; % Expected number of pathces (not exactly value)

% Point charge
pointchrg = -1; % q=1e
xp = 1.5*Radius;
yp = 0;
zp = 0;

% Note: Theta and Phi are used oppositely compared to the 1-page paper
a = (4*pi*Radius^2)/N; % Average area per patch
d = sqrt(a); % Length of patch
M_Phi = round(pi*Radius/d);
d_Phi = pi*Radius/M_Phi;
d_Theta = a/d_Phi;

% The center position of each patch
xc = linspace(0,0,N)';
yc = linspace(0,0,N)';
zc = linspace(0,0,N)';
% Normal unit vector
NormVec_x = linspace(0,0,N)';
NormVec_y = linspace(0,0,N)';
NormVec_z = linspace(0,0,N)';
% Patch area
DeltaArea = linspace(0,0,N)';

count = 0;
for i = 0:(M_Phi-1)
    
    Phi = pi*(i+0.5)/M_Phi;
    M_Theta = round(2*pi*Radius*sin(Phi)/d_Theta);
    
    for j = 0:(M_Theta-1)
        count = count + 1;
        
        Theta = 2*pi*j/M_Theta;
        xc(count)=Radius*sin(Phi)*cos(Theta);
        yc(count)=Radius*sin(Phi)*sin(Theta);
        zc(count)=Radius*cos(Phi);
    
        NormVec_x(count) = sin(Phi)*cos(Theta);
        NormVec_z(count) = sin(Phi)*sin(Theta);
        NormVec_y(count) = cos(Phi);

    end
end

% Approximated patch area
NN = length(xc); % Number of generated patches
DeltaArea = (4*pi*Radius^2)/NN * linspace(1,1,NN)'; % Average area (should be improved)

% This is used to check the location of patches
% figure(1);
% size = 10;
% scatter3(xc(:),yc(:),zc(:),size,'filled');
% hold on;
% 
% scatter3(xp,yp,zp,size,'filled');
% axis equal;
TotalArea = sum(DeltaArea,'all');

% Electrical field by free charge
Ex_f = linspace(0,0,NN)';
Ey_f = linspace(0,0,NN)';
Ez_f = linspace(0,0,NN)';
dotEn_f = linspace(0,0,NN)'; % Dot product: E*n

for i = 1:NN
    
    % I removed the contribution by free surface charge here
    % Sigma_f = 0
    
    % Contribution by point charge
    dx = xc(i)-xp;
    dy = yc(i)-yp;
    dz = zc(i)-zp;

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);

    dEx_p = dx/(dr*dr2)*pointchrg/4/pi/epsilon;
    dEy_p = dy/(dr*dr2)*pointchrg/4/pi/epsilon;
    dEz_p = dz/(dr*dr2)*pointchrg/4/pi/epsilon;

    Ex_f(i) = Ex_f(i) + dEx_p;
    Ey_f(i) = Ey_f(i) + dEy_p;
    Ez_f(i) = Ez_f(i) + dEz_p;
    
    % Compute the dot product term
    dotEn_f(i) = Ex_f(i)*NormVec_x(i) + Ey_f(i)*NormVec_y(i)+...
            Ez_f(i)*NormVec_z(i);
end

% Iterate to find bound charge density
kappa_ave = 0.5*(kappa_air+kappa_p);
kappa_diff = kappa_air-kappa_p; % kappa_out-kappa_in

b = (1-kappa_ave)*sigma_f - kappa_diff*epsilon*dotEn_f;

% Construct L_i,j
L = zeros(NN);
% Electrical field by bound charge
dotEn_b = linspace(0,0,NN)'; % Dot product: E*n
    
    for i = 1:NN
    for j = 1:NN
        
        if (i==j)
            
            L(i,j)=kappa_ave; % the patch itself
            
        else %% if i~=j
        
        dx = xc(i)-xc(j);
        dy = yc(i)-yc(j);
        dz = zc(i)-zc(j);
        
        dr2 = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr2);
        
        dEx_b = dx/(dr*dr2)*DeltaArea(j)/4/pi;
        dEy_b = dy/(dr*dr2)*DeltaArea(j)/4/pi;
        dEz_b = dz/(dr*dr2)*DeltaArea(j)/4/pi;
        
        dotEn_b(i) = dEx_b*NormVec_x(i) + dEy_b*NormVec_y(i)+...
            dEz_b*NormVec_z(i);
        L(i,j) = dotEn_b(i);
        end % if (i==j)
    end % for j=1:NN
    end % for i=1:NN

sigma_b = L\b;
error = norm(L*sigma_b-b);
fprintf('Error = %f\n',error);

sigma = sigma_f + sigma_b;

figure(2);
size = 10;
scatter3(xc(:),yc(:),zc(:),size,sigma(:),'filled');
hold on;
scatter3(xp,yp,zp,size,'r','filled');
axis equal;

% Calculate potential at point charge location
potential = linspace(0,0,NN)';

for i = 1:NN
    %Potential by surface charge
    for j = 1:NN
    
    if (i==j)
        continue;
    end
        
    dx = xc(i) - xc(j);
    dy = yc(i) - yc(j);
    dz = zc(i) - zc(j);

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);
    
    potential(i) = potential(i) + sigma(j)*DeltaArea(j)/dr/4/pi/epsilon;
    
    end
    
    %Potential by point charge
    dx = xc(i)-xp;
    dy = yc(i)-yp;
    dz = zc(i)-zp;

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);
    
    potential(i) = potential(i) + pointchrg/dr/4/pi/epsilon;
    
end

distance = xp/Radius - 1;
E_Patch = (sigma.*potential).*DeltaArea;
Energy = 0.5*sum(E_Patch);
NormEnergy = Energy/(pointchrg^2)/epsilon/Radius^2;
fprintf('Surface distance = %f\n',distance);
fprintf('Normalized Energy = %f\n',Energy);
% total_sigma = sum(sigma);
% fprintf('Total sigma = %f\n',total_sigma);
% total_potential = sum(potential);
% fprintf('Total potential = %f\n',total_potential);
% 
% figure(3);
% scatter(potential,sigma);