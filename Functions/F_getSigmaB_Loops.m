function [sigma_b,b] = F_getSigmaB_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj)
% PROVIDES VECTOR OF BOUND CHARGE SURFACE DENSITIES FOR EACH PATCH
%{   
    Given:
    R..................... sphere radius
    x,y,z................. locations of the patches
    nVect................. matrix of normal vectors for each patch (nVect(i) = nVx(i), nVy(i), nVz(i))
    x_pc,y_pc,z_pc........ location of the point charge
    pcharge............... charge of point charge
    sigma_f............... vector of patch free charge surface densities 
    k_air................. dielectric constant of the air (or surrounding medium)
    k_obj................. dielectric constant of the sphere material

    Returns;
    sigma_b............... vector of patch bound charge surface densities
    b..................... vector of known charge (A*sigma_b = b). See Barrios and Luijten 2014, Journ. of Chem Phys. 
    
%}

Npatches = length(x);
dA = 4*pi*(R^2)/Npatches;

A1 = zeros(Npatches);
B1 = zeros(Npatches);

k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);

% Patch contributions
for i = 1:Npatches
for j = 1:Npatches
    
    dx = x(i)-x(j); dy = y(i)-y(j); dz = z(i)-z(j);
    dr = sqrt(dx^2 + dy^2 + dz^2);
    
    if i==j
        A1(i,j) = k_bar;
        B1(i,j) = 1-k_bar;
    else
        A1(i,j) = k_delta*dA/4/pi*(dx*nVect(i,1) + dy*nVect(i,2) + dz*nVect(i,3))/(dr^3);
        B1(i,j) = -k_delta*dA/4/pi*(dx*nVect(i,1) + dy*nVect(i,2) + dz*nVect(i,3))/(dr^3);
    end
end
end

% Point charge contribution
b_pc = zeros(Npatches,1);
for i = 1:Npatches
    dx_pc = x(i)-x_pc; dy_pc = y(i)-y_pc; dz_pc = z(i)-z_pc;
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    
    b_pc(i) = (dx_pc*nVect(i,1) + dy_pc*nVect(i,2) + dz_pc*nVect(i,3))/(dr_pc^3);
end

b = B1*sigma_f - k_delta*pcharge/4/pi*(b_pc);
sigma_b = A1\b;
end

