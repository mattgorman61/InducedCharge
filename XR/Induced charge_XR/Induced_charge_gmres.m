clear;
clc;

% ========== Definition ==========
sigma_f = 0; % Surface free charge density
kappa_air = 1; % Dielectric constance of air
kappa_p = 0.1; % Dielectric constance of target
epsilon = 1; %Unit: F/m

% ========== Geometry ==========
% Geometry: sphere + point charge
Radius = 1; % Sphere radius: 100 micron
N = 376; % Expected number of pathces (not exactly value)

% Point charge
pointchrg = -1;
xp = 1.5*Radius;
yp = 0;
zp = 0;

% Generate patch geometry
[xc,yc,zc,NormVec,DeltaArea,NN]=F_GeometrySphere(Radius,N);

% This is used to check the location of patches
% figure(1);
% size = 10;
% scatter3(xc(:),yc(:),zc(:),size,'filled');
% axis equal;
% 
% TotalArea = sum(DeltaArea,'all');

rdiff = zeros(NN,NN,3);
rdiff(:,:,1) = xc - xc';
rdiff(:,:,2) = yc - yc';
rdiff(:,:,3) = zc - zc';

rdist1 = pdist([xc,yc,zc]);
rdist = squareform(rdist1);

I = rdiff ./(rdist.^3)/4/pi; % Diagonal components become NaN
for i = 1:NN
    I(i,i,:) = [0,0,0];
end

% Check_I_1=I(:,:,1);
% Check_I_2=I(:,:,2);
% Check_I_3=I(:,:,3);

L = zeros(NN,NN);
for i = 1:NN
    L(i,:) = I(i,:,1)*NormVec(i,1)+I(i,:,2)*NormVec(i,2)+I(i,:,3)*NormVec(i,3);
end

k_ave = 0.5*(kappa_air+kappa_p);
k_diff = kappa_air-kappa_p; %kappa_out - kappa_in

dSj = DeltaArea.*linspace(1,1,NN)';% Assume patch area is the same,can be updated if necessary
A = k_ave*eye(NN)+k_diff*L*dSj(1); % Left matrix
% Check_A = k_diff*L*dSj(1);

sig_f = sigma_f * linspace(1,1,NN)';
B = (1-k_ave)*eye(NN)-k_diff*L*dSj(1); % Geometry matrix for free charge

% Compute contribution by point charge
xdist = xc - xp;
ydist = yc - yp;
zdist = zc - zp;

pointchrgdist = sqrt(xdist.^2 + ydist.^2 + zdist.^2);
Ex_ext = pointchrg*xdist./(pointchrgdist.^3)/4/pi;
Ey_ext = pointchrg*ydist./(pointchrgdist.^3)/4/pi;
Ez_ext = pointchrg*zdist./(pointchrgdist.^3)/4/pi;

E_proj=linspace(0,0,NN)';
for i = 1:NN
E_proj(i)=Ex_ext(i)*NormVec(i,1)+Ey_ext(i)*NormVec(i,2)+Ez_ext(i)*NormVec(i,3);
end

b = B*sig_f - k_diff*E_proj;

sig_b = gmres(A,b,10,1E-4); % Tolerance: 1E-4

sigma = sig_f + sig_b;

figure(2);
size = 10;
scatter3(xc(:),yc(:),zc(:),size,sigma(:),'filled');
hold on;
scatter3(xp,yp,zp,size,'r','filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gca,'LineWidth',1.5);

% Compute the energy
Phi_Patch2Point = DeltaArea.*sigma./(pointchrgdist)/4/pi/epsilon;
Phi_Point_Total = sum(Phi_Patch2Point);
U_Point = 0.5 * Phi_Point_Total * pointchrg;
U0 = pointchrg^2/(epsilon*kappa_air*Radius);
U_Norm = U_Point/U0;
fprintf('Computed normalized energy is: %f\n',U_Norm);

net_chrg = sum(sigma);
fprintf('Net charge is: %f\n',net_chrg);


% % Check results
% potential = linspace(0,0,NN)';
% 
% for i = 1:NN
%     %Potential by surface charge
%     for j = 1:NN
%     
%     if (i==j)
%         continue;
%     end
%         
%     dx = xc(i) - xc(j);
%     dy = yc(i) - yc(j);
%     dz = zc(i) - zc(j);
% 
%     dr2 = dx*dx + dy*dy + dz*dz;
%     dr = sqrt(dr2);
%     
%     potential(i) = potential(i) + sigma(j)*DeltaArea(j)/dr/4/pi/epsilon;
%     
%     end
%     
%     %Potential by point charge
%     dx = xc(i)-xp;
%     dy = yc(i)-yp;
%     dz = zc(i)-zp;
% 
%     dr2 = dx*dx + dy*dy + dz*dz;
%     dr = sqrt(dr2);
%     
%     potential(i) = potential(i) + pointchrg/dr/4/pi/epsilon;
%     
% end

% % Compute theoratical value: From Jones, 1995.
%     N_th = 100;
%     N_term = 1;
%     angle_th = linspace(0,pi,N_th)';
%     % Coordinates of points (x,y,0)
%     x_th = Radius*cos(angle_th)';
%     y_th = Radius*sin(angle_th)';
%     z_th = linspace(0,0,N_th)';
%     potential_th = linspace(0,0,N_th)';
%     for i = 1:N_th
%         dx = xp-x_th(i);
%         dy = yp-y_th(i);
%         dz = zp-z_th(i);
%         dr = sqrt(dx*dx+dy*dy+dz*dz);
%         % Contribution by point charge
%         potential_th(i) = potential_th(i)+pointchrg/4/pi/epsilon/dr;
%         % Contribution by higher-order terms
%         k_ratio = kappa_p/kappa_air;
%         for j = 0:N_term
%             % Compute coefficient A
%             A = -(pointchrg/4/pi/epsilon)*(j*(k_ratio-1)/(j*(k_ratio+1)+1))...
%                 *(Radius^(2*j+1))/((xp)^(j+1));
%             % Compute Ledengre polynomials
%             P = 0;
%             x = cos(angle_th(i));
%             for k = 0:j
%                 P = P + (nchoosek(j,k)^2)*((x-1)^(j-k))*((x+1)^(k));
%             end
%             P = P/(2^j);
%             potential_th(i) = potential_th(i) + A*P/(Radius^(j+1));
%         end
%     end
%     
%     % Plot potential profile of z=0 for our calculation
%     angle = acos(xc(2936:3004)/Radius);% 2936:3004 is only for N = 6000
%     figure(3);
%     p1=plot(angle,potential(2936:3004),'x','Linewidth',1.5);
%     hold on;
%     p2=plot(angle_th,potential_th,'Linewidth',1.5);
%     set(gca,'LineWidth',1.5);
%     xlabel('\theta','Fontsize',14);
%     ylabel('Potential','Fontsize',14);
%     
%     h=legend([p1,p2],'Our calculation','Theory','Fontsize',14);
%     set(h,'box',' off');