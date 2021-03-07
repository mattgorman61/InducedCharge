close all;
clear all;
clc;

% Definition
sigma_f = 0; % Surface free charge density
kappa_air = 1; % Dielectric constance of air
kappa_p = 3.6; % Dielectric constance of target
itmax = 5; % Maximum iteration

% Geometry: sphere + point charge
Theta = 2*pi;% Angle of Longitude
Phi = pi;% Angle of Latitude
Radius = 5;

NTheta = 15; % Longitude
NPhi = 15; % Latitude
DeltaTheta = Theta/NTheta;
DeltaPhi = Phi/NPhi;

pointchrg = -1;
xp = 1.1*Radius;
yp = 0;
zp = 0;

% The center position of each patch
xc = zeros(NTheta,NPhi);
yc = zeros(NTheta,NPhi);
zc = zeros(NTheta,NPhi);
% Normal unit vector
NormVec_x = zeros(NTheta,NPhi);
NormVec_y = zeros(NTheta,NPhi);
NormVec_z = zeros(NTheta,NPhi);
% Patch area
DeltaArea = zeros(NTheta,NPhi);

for i = 1:NTheta
for j = 1:NPhi

    xc(i,j)=Radius*sin(DeltaPhi*(j-0.5))*cos(DeltaTheta*(i-0.5));
    yc(i,j)=Radius*sin(DeltaPhi*(j-0.5))*sin(DeltaTheta*(i-0.5));
    zc(i,j)=Radius*cos(DeltaPhi*(j-0.5));
    
    NormVec_x(i,j) = sin(DeltaPhi*(i-0.5))*cos(DeltaTheta*(j-0.5));
    NormVec_z(i,j) = sin(DeltaPhi*(i-0.5))*sin(DeltaTheta*(j-0.5));
    NormVec_y(i,j) = cos(DeltaPhi*(j-0.5));
    
    DeltaArea(i,j) = abs(Radius*Radius*sin(DeltaTheta*(i-0.5))*DeltaTheta*DeltaPhi);
    
end
end

% figure(1);
% size = 10;
% scatter3(xc(:,1),yc(:,1),zc(:,1),size,'filled');
% hold on;
% for i = 2:NPhi
%     scatter3(xc(:,i),yc(:,i),zc(:,i),size,'filled');
% end
% scatter3(xp,yp,zp,size,'filled');
% axis equal;
TotalArea = sum(DeltaArea,'all');

% Electrical field by free charge
Ex_f = zeros(NTheta,NPhi);
Ey_f = zeros(NTheta,NPhi);
Ez_f = zeros(NTheta,NPhi);
dotEn_f = zeros(NTheta,NPhi); % Dot product: E*n

for i = 1:NTheta
for j = 1:NPhi
    % Contribution of cylinder patches
    for m = 1:NTheta
    for n = 1:NPhi
        
        if ((m==i)&&(n==j))
            continue; % skip the patch itself
        end
        
        dx = xc(i,j)-xc(m,n);
        dy = yc(i,j)-yc(m,n);
        dz = zc(i,j)-zc(m,n);
        
        dr2 = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr2);
        
        dEx_f = dx/(dr*dr2)*DeltaArea(i,j)*sigma_f;
        dEy_f = dy/(dr*dr2)*DeltaArea(i,j)*sigma_f;
        dEz_f = dz/(dr*dr2)*DeltaArea(i,j)*sigma_f;
        
        Ex_f(i,j) = Ex_f(i,j) + dEx_f;
        Ey_f(i,j) = Ey_f(i,j) + dEy_f;
        Ez_f(i,j) = Ez_f(i,j) + dEz_f;
        
        checkmax = max(max(dEx_f,dEy_f),dEz_f);
        checkmin = min(min(dEx_f,dEy_f),dEz_f);
        
        if ((checkmax>10)||(checkmin<(-10)))
            fprintf('i=%d, j=%d, m=%d, n=%d\n',i,j,m,n);
            fprintf('dEx_f=%f,dEy_f=%f,dEz_f=%f\n',dEx_f,dEy_f,dEz_f);
            fprintf('xc(%d,%d)=%f,yc(%d,%d)=%f,zc(%d,%d)=%f,xc(%d,%d)=%f,yc(%d,%d)=%f,zc(%d,%d)=%f\n',...
                i,j,xc(i,j),i,j,yc(i,j),i,j,zc(i,j),m,n,xc(m,n),m,n,yc(m,n),m,n,zc(m,n));
%             pause(10);
        end
        
    end
    end
    
    % Contribution point charge
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
    Ex_b = zeros(NTheta,NPhi);
    Ey_b = zeros(NTheta,NPhi);
    Ez_b = zeros(NTheta,NPhi);
    dotEn_b = zeros(NTheta,NPhi); % Dot product: E*n
    
    for i = 1:NTheta
    for j = 1:NPhi
    for m = 1:NTheta
    for n = 1:NPhi
        
        if ((m==i)&&(n==j))
            continue; % skip the patch itself
        end
        
        dx = xc(i,j)-xc(m,n);
        dy = yc(i,j)-yc(m,n);
        dz = zc(i,j)-zc(m,n);
        
        dr2 = dx*dx + dy*dy + dz*dz;
        dr = sqrt(dr2);
        
        dEx_b = dx/(dr*dr2)*DeltaArea(i,j)*sigma_b(m,n);
        dEy_b = dy/(dr*dr2)*DeltaArea(i,j)*sigma_b(m,n);
        dEz_b = dz/(dr*dr2)*DeltaArea(i,j)*sigma_b(m,n);
        
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

figure(2);
size = 10;
scatter3(xc(:,1),yc(:,1),zc(:,1),size,sigma(:,1),'filled');
hold on;
for i = 2:NPhi
    scatter3(xc(:,i),yc(:,i),zc(:,i),size,sigma(:,i),'filled');
end
scatter3(xp,yp,zp,size,'r','filled');
axis equal;

% Calculate potential at point charge location
potential = 0;
for i = 1:NTheta
for j = 1:NPhi
    
    % Contribution to point charge
    dx = xp - xc(i,j);
    dy = yp - yc(i,j);
    dz = zp - zc(i,j);

    dr2 = dx*dx + dy*dy + dz*dz;
    dr = sqrt(dr2);
    
    potential = potential + sigma(i,j)*DeltaArea(i,j)/dr;

end
end

distance = xp - Radius;
Energy = 0.5*potential*pointchrg;
fprintf('Surface distance = %f\n',distance);
fprintf('Normalized Energy = %f\n',Energy);
