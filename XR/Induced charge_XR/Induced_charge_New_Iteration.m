clear;
clc;

% ========== Definition ==========
sigma_f = 0; % Surface free charge density
kappa_air = 1; % Dielectric constance of air
kappa_p = 2.5; % Dielectric constance of target
epsilon = 1; %Unit: F/m

% ========== Geometry ==========
% Geometry: sphere + point charge
Radius = 1; % Sphere radius: 100 micron
N = 2000; % Expected number of pathces (not exactly value)

% Point charge
pointchrg = -1;
xp = 1.5*Radius;
yp = 0;
zp = 0;

% Generate patch geometry
[xc,yc,zc,NormVec,DeltaArea,NN]=F_GeometrySphere(Radius,N);
% quiver3(xc,yc,zc,NormVec(:,1),NormVec(:,2),NormVec(:,3),1);
% axis equal;

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

rdiff_Check_1 = xc - xc';
rdiff_Check_2 = yc - yc';
rdiff_Check_3 = zc - zc';

rdist1 = pdist([xc,yc,zc]);
rdist = squareform(rdist1);

I = rdiff ./(rdist.^3)/4/pi; % Diagonal components become NaN
for i = 1:NN
    I(i,i,:) = [0,0,0];
end

I_Check_1 = I(:,:,1);
I_Check_2 = I(:,:,2);
I_Check_3 = I(:,:,3);

L = zeros(NN,NN);
for i = 1:NN
    for j = 1:NN
    L(i,j) = I(i,j,1)*NormVec(i,1)+I(i,j,2)*NormVec(i,2)+I(i,j,3)*NormVec(i,3);
    end
end

k_ave = 0.5*(kappa_air+kappa_p);
k_diff = kappa_air-kappa_p; %kappa_out - kappa_in

dSj = DeltaArea.*linspace(1,1,NN)';% Assume patch area is the same,can be updated if necessary
A = k_ave*eye(NN)+k_diff*L*dSj(1); % Left matrix


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

% Use Richardson iteration to solve the equation
sig_b_old = linspace(0,0,NN)';
RichCount = 1;
norm_b = norm(b,2);
while(1)
        
    % Iteration: I use sig_b and sig_b_old for clarity, actually not
    % necessary
    sig_b = sig_b_old + 1/k_ave*(b-A*sig_b_old);
    sig_b_old = sig_b;
    
    % Compute the error
    error = b - A*sig_b;
    norm_error = norm(error,2);
    
    if (norm_error*10000 < norm_b)
        fprintf('Richardson Iteration converges at: %d\n',RichCount);
        fprintf('2-Norm of error is: %f\n',norm_error);
        break;
    end
    RichCount = RichCount + 1;
end

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

% Extract potential info on circle yc=0
l_PotentialCircle = (yc==0);
PotentialCircle = linspace(0,0,NN)';

for i = 1:NN
    if (l_PotentialCircle(i)==1)
        for j = 1:NN
            if (i==j)
                continue; % Jump itself
            end
            PotentialCircle(i) = PotentialCircle(i) + sigma(j)*DeltaArea(j)/4/pi/epsilon/rdist(i,j);
        end
        PotentialCircle(i) = PotentialCircle(i) + pointchrg/4/pi/epsilon/pointchrgdist(i);
    end
end

NCircle = sum(l_PotentialCircle);
PotentialCirclePlot = linspace(0,0,NCircle/2)';
AnglePlot = linspace(0,0,NCircle/2)';
CircleCount = 0;
angle = acos(xc/Radius);

for i = 1:NN
    if((l_PotentialCircle(i)==1)&&(zc(i)>=0))
        CircleCount = CircleCount + 1;
        AnglePlot(CircleCount) = angle(i);
        PotentialCirclePlot(CircleCount) = PotentialCircle(i);
    end
end

% Compute theoratical value: From Jones, 1995.
% This theoratical results seem to fluctuate at large kappa_p/kappa_m

N_th = 100;
N_term = 20;
angle_th = linspace(0,pi,N_th)';
% Coordinates of points (x,y,0)
x_th = Radius*cos(angle_th)';
y_th = Radius*sin(angle_th)';
z_th = linspace(0,0,N_th)';
potential_th = linspace(0,0,N_th)';

for i = 1:N_th
    dx = xp-x_th(i);
    dy = yp-y_th(i);
    dz = zp-z_th(i);
    dr = sqrt(dx*dx+dy*dy+dz*dz);
    % Contribution by point charge
    potential_th(i) = potential_th(i)+pointchrg/4/pi/epsilon/dr;
    % Contribution by higher-order terms
    k_ratio = kappa_p/kappa_air;
    for j = 0:N_term
        % Compute coefficient An
        An = -(pointchrg/4/pi/epsilon)*(j*(k_ratio-1)/(j*(k_ratio+1)+1))...
            *(Radius^(2*j+1))/((xp)^(j+1));
        % Compute Ledengre polynomials
        P = 0;
        x = cos(angle_th(i));
        for k = 0:j
            P = P + (nchoosek(j,k)^2)*((x-1)^(j-k))*((x+1)^(k));
        end
        P = P/(2^j);
        potential_th(i) = potential_th(i) + An*P/(Radius^(j+1));
    end
end

% Plot potential profile comparison
Potential_0 = abs(pointchrg)/4/pi/epsilon/Radius; % For normalization
potential_th = potential_th/Potential_0;
PotentialCirclePlot = PotentialCirclePlot/Potential_0;

figure(3);
p2=plot(angle_th,potential_th,'Linewidth',1.5);
hold on;
p1 = plot(AnglePlot,PotentialCirclePlot,'x');

set(gca,'LineWidth',1.5);
xlabel('\theta','Fontsize',14);
ylabel('Potential','Fontsize',14);

h=legend([p1,p2],'Our calculation','Theory','Fontsize',14);
set(h,'box',' off');
