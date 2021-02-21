clear;
clc;

a = 4;
b = 4;
c = 1;

kappa_p = 2.5;
kappa_air = 1;

Ex = 1;
Ey = 1;
Ez = 1;

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
factor1 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Lx);
factor2 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Ly);
factor3 = kappa_air*(k_ratio-1)/(1+(k_ratio-1)*Lz);

px_eff = 4*pi*a*b*c/3*factor1*Ex;
py_eff = 4*pi*a*b*c/3*factor2*Ey;
pz_eff = 4*pi*a*b*c/3*factor3*Ez;

fprintf('px_eff=%f, py_eff=%f, pz_eff=%f\n',px_eff,py_eff,pz_eff);

% Compute theoratical torque
coeff = 4*pi*a*b*c*((kappa_p-kappa_air)^2)/3/kappa_air;

Tx = coeff*(Lz-Ly)*Ey*Ez/(1+(k_ratio-1)*Ly)/(1+(k_ratio-1)*Lz);
Ty = coeff*(Lx-Lz)*Ex*Ez/(1+(k_ratio-1)*Lx)/(1+(k_ratio-1)*Lz);
Tz = coeff*(Ly-Lx)*Ey*Ex/(1+(k_ratio-1)*Ly)/(1+(k_ratio-1)*Lx);

fprintf('Tx=%f, Ty=%f, Tz=%f\n',Tx,Ty,Tz);

