clear;
clc;

kappa_p = 2.5;
kappa_air = 1;
k_ratio = kappa_p/kappa_air;

Radius = 1;
E = 1;

p_eff = 4*pi*kappa_air*(k_ratio-1)/(k_ratio+2)*(Radius^3)*E;

fprintf('Induced dipole of the sphere = %f\n',p_eff);

