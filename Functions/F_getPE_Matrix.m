function [U_pCharge_Norm,U_tot_Norm,netCharge] = F_getPE_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0)
% PROVIDES NORMALIZED POTENTIAL ENERGY OF THE POINT CHARGE SPHERE SYSTEM
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

U_tot_patches = 0;
U_tot_pcharge = 0;

k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);

for i = 1:Npatches
% Patch-patch contribution
for j = 1:Npatches
    if (i~=j)
    dx = x(i)-x(j); dy = y(i)-y(j); dz = z(i)-z(j);
    dr = sqrt(dx^2 + dy^2 + dz^2);
    U_tot_patches = U_tot_patches + (sigma(i)*(dA^2))/(dr);
    end
end
% Pcharge-patch contribution
    dx_pc = x(i)-x_pc;  dy_pc = y(i)-y_pc; dz_pc = z(i)-z_pc;
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    U_tot_patches = U_tot_patches + sigma(i)*dA*pcharge/(dr_pc);
end

% Patch-pcharge contribution
for i = 1:Npatches
    dx_pc = x_pc-x(i);  dy_pc = y_pc-y(i); dz_pc = z_pc-z(i);
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    U_tot_pcharge = U_tot_pcharge + sigma(i)*dA*pcharge/(dr_pc);
end

U_tot = 0.5*(1/4/pi/epsilon_0)*(U_tot_patches + U_tot_pcharge);
U0 = (pcharge^2)/epsilon_0/k_air/R;
U_tot_Norm = U_tot/U0;

U_pCharge_Norm =  (0.5/4/pi/epsilon_0) * U_tot_pcharge / U0;
netCharge = sum(sigma);

%{
fprintf('RESULTS:\n')
fprintf('k_tilda: %.5f\n', k_tilda);
fprintf('Surface Distance: %f\n', surfDist);
fprintf('Normalized PE: %f\n', U_tot_Norm);
fprintf('Net Charge: %f\n\n\n', netCharge);
%}

end

