clear all;
close all;
clc;

% Definition
sigma_f = 1; % Surface free charge density
kappa_air = 1; % Dielectric constance of air
kappa_p = 3.6; % Dielectric constance of target
itmax = 3; % Maximum iteration

% Geometry: cylinder + negative point charge
Width = 5;
Radius = 1;

NCircle = 100;
NWidth = 100;
DeltaTheta = 2*pi/NCircle;
DeltaWidth = Width/NWidth;
DeltaArea = Radius*sin(DeltaTheta/2)*DeltaWidth;

pointchrg = 1;
xp = 2*Radius;
yp = 0.5*Width;
zp = 0.5*Radius;

% The center position of each patch
xc = zeros(NCircle,NWidth);
yc = zeros(NCircle,NWidth);
zc = zeros(NCircle,NWidth);
% Normal unit vector
NormVec_x = zeros(NCircle,NWidth);
NormVec_y = zeros(NCircle,NWidth);
NormVec_z = zeros(NCircle,NWidth);

for i = 1:NCircle
for j = 1:NWidth

    xc(i,j)=Radius*cos(DeltaTheta*(i-0.5));
    yc(i,j)=DeltaWidth*(j-0.5);
    zc(i,j)=Radius*sin(DeltaTheta*(i-0.5));
    
    NormVec_x(i,j) = cos(DeltaTheta*(i-0.5));
    NormVec_z(i,j) = sin(DeltaTheta*(i-0.5));
    NormVec_y(i,j) = 0; % Normal vector lie in xOz plane
    
end
end

% Electrical field by free charge
Ex_f = zeros(NCircle,NWidth);
Ey_f = zeros(NCircle,NWidth);
Ez_f = zeros(NCircle,NWidth);
dotEn_f = zeros(NCircle,NWidth); % Dot product: E*n

for i = 1:NCircle
for j = 1:NWidth
    % Contribution of cylinder patches
    for m = 1:NCircle
    for n = 1:NWidth
        
        if ((m==i)&&(n==j))
            continue; % skip the patch itself
        end
        
        dx = xc(i,j)-xc(m,n);
        dy = yc(i,j)-yc(m,n);
        dz = zc(i,j)-zc(m,n);
        
        dr2 = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr2);
        
        dEx_f = dx/(dr*dr2)*DeltaArea*sigma_f;
        dEy_f = dy/(dr*dr2)*DeltaArea*sigma_f;
        dEz_f = dz/(dr*dr2)*DeltaArea*sigma_f;
        
        Ex_f(i,j) = Ex_f(i,j) + dEx_f;
        Ey_f(i,j) = Ey_f(i,j) + dEy_f;
        Ez_f(i,j) = Ez_f(i,j) + dEz_f;
        
    end
    end
    
    % Contribution of negative point charge
    dx = xc(i,j)-xp;
    dy = yc(i,j)-yp;
    dz = zc(i,j)-zp;

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);

    dEx_p = dx/(dr*dr2)*pointchrg;
    dEy_p = dy/(dr*dr2)*pointchrg;
    dEz_p = dz/(dr*dr2)*pointchrg;

    Ex_f(i,j) = Ex_f(i,j) + dEx_p;
    Ey_f(i,j) = Ey_f(i,j) + dEy_p;
    Ez_f(i,j) = Ez_f(i,j) + dEz_p;
    
    % Compute the dot product term
    dotEn_f(i,j) = Ex_f(i,j)*NormVec_x(i,j) + Ey_f(i,j)*NormVec_y(i,j)+...
            Ez_f(i,j)*NormVec_z(i,j);
end
end

% Iterate to find bound charge density
kappa_ave = 0.5*(kappa_air+kappa_p);
kappa_diff = kappa_air-kappa_p; % kappa_out-kappa_in

b = (1-kappa_ave)*sigma_f - kappa_diff/4/pi*dotEn_f;
sigma_b = b/kappa_ave;

for it = 1:itmax
    % Electrical field by bound charge
    Ex_b = zeros(NCircle,NWidth);
    Ey_b = zeros(NCircle,NWidth);
    Ez_b = zeros(NCircle,NWidth);
    dotEn_b = zeros(NCircle,NWidth); % Dot product: E*n
    
    for i = 1:NCircle
    for j = 1:NWidth
    for m = 1:NCircle
    for n = 1:NWidth
        
        if ((m==i)&&(n==j))
            continue; % skip the patch itself
        end
        
        dx = xc(i,j)-xc(m,n);
        dy = yc(i,j)-yc(m,n);
        dz = zc(i,j)-zc(m,n);
        
        dr2 = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr2);
        
        dEx_b = dx/(dr*dr2)*DeltaArea*sigma_b(m,n);
        dEy_b = dy/(dr*dr2)*DeltaArea*sigma_b(m,n);
        dEz_b = dz/(dr*dr2)*DeltaArea*sigma_b(m,n);
        
        Ex_b(i,j) = Ex_b(i,j) + dEx_b;
        Ey_b(i,j) = Ey_b(i,j) + dEy_b;
        Ez_b(i,j) = Ez_b(i,j) + dEz_b;
        
        dotEn_b(i,j) = Ex_b(i,j)*NormVec_x(i,j) + Ey_b(i,j)*NormVec_y(i,j)+...
            Ez_b(i,j)*NormVec_z(i,j);
    end
    end
    end
    end
    %End of compute E induced by bound charge
    sigma_b = (b - kappa_diff/4/pi*dotEn_b)/kappa_ave;
end

sigma = sigma_f + sigma_b;

figure(1);
size = 10;
scatter3(xc(:,1),yc(:,1),zc(:,1),size,sigma(:,1),'filled');
hold on;
for i = 2:NWidth
    scatter3(xc(:,i),yc(:,i),zc(:,i),size,sigma(:,i),'filled');
end
scatter3(xp,yp,zp,size,'filled');
axis equal;
colorbar;

% Potential at point charge location
potential = 0;
for i = 1:NCircle
for j = 1:NWidth
        
    % Contribution of negative point charge
    dx = xp-xc(i,j);
    dy = yp-yc(i,j);
    dz = zp-zc(i,j);

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);

    dEx_p = dx/(dr*dr2)*pointchrg;
    dEy_p = dy/(dr*dr2)*pointchrg;
    dEz_p = dz/(dr*dr2)*pointchrg;

    Ex_f(i,j) = Ex_f(i,j) + dEx_p;
    Ey_f(i,j) = Ey_f(i,j) + dEy_p;
    Ez_f(i,j) = Ez_f(i,j) + dEz_p;
    
end
end
