clear;
clc;

a = 4E-4; % m
b = 1E-4; % m
c = 1E-4; % m

kappa_p = 2.5;
kappa_air = 1;
epsilon = 8.854E-12; % F/m

Ex = 1E5; %F/m
Ey = 1E5; %F/m
Ez = 1E5; %F/m

k_ratio = kappa_p/kappa_air;

% Elliptical integral
syms s
f1 = 1/(s+a^2)^(1.5)/(s+b^2)^(0.5)/(s+c^2)^(0.5);
F1 = int(f1,s,0,+inf);
Lx = double(F1);

f2 = 1/(s+a^2)^(0.5)/(s+b^2)^(1.5)/(s+c^2)^(0.5);
F2 = int(f2,s,0,+inf);
Ly = double(F2);

f3 = 1/(s+a^2)^(0.5)/(s+b^2)^(0.5)/(s+c^2)^(1.5);
F3 = int(f3,s,0,+inf);
Lz = double(F3);

Lx = Lx*a*b*c/2;
Ly = Ly*a*b*c/2;
Lz = Lz*a*b*c/2;

fprintf('Sum of ellip int = %f\n',Lx+Ly+Lz);

% Compute effective dipole
% factor1 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Lx);
% factor2 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Ly);
% factor3 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Lz);
% 
% px_eff = 4*pi**a*b*c/3*factor1*Ex;
% py_eff = 4*pi*a*b*c/3*factor2*Ey;
% pz_eff = 4*pi*a*b*c/3*factor3*Ez;

% fprintf('px_eff=%f, py_eff=%f, pz_eff=%f\n',px_eff,py_eff,pz_eff);

% Compute theoratical torque
coeff = 4*pi*epsilon*a*b*c*((kappa_p-kappa_air)^2)/3/kappa_air;

Tx = coeff*(Lz-Ly)*Ey*Ez/(1+(k_ratio-1)*Ly)/(1+(k_ratio-1)*Lz); %Unit: [N*m]
Ty = coeff*(Lx-Lz)*Ex*Ez/(1+(k_ratio-1)*Lx)/(1+(k_ratio-1)*Lz); %Unit: [N*m]
Tz = coeff*(Ly-Lx)*Ey*Ex/(1+(k_ratio-1)*Ly)/(1+(k_ratio-1)*Lx); %Unit: [N*m]

Tx = Tx * 1E15; % Unit: [nN*um]
Ty = Ty * 1E15; % Unit: [nN*um]
Tz = Tz * 1E15; % Unit: [nN*um]
fprintf('Tx=%f nN*um\n',Tx);
fprintf('Ty=%f nN*um\n',Ty);
fprintf('Tz=%f nN*um\n',Tz);

