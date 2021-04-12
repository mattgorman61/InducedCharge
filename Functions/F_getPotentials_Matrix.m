function [phi,phi_0,phi_norm,theta] = F_getPotentials_Matrix(R,x,y,z,dA,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0)
% PROVIDES VECTOR OF FORCES ACTING ON EACH PATCH OF THE SPHERES
%{   
    Given:
    R..................... sphere radius
    x,y,z................. locations of the patches
    dA ................... vector of patch areas, constant for spheres.
    nVect................. matrix of normal vectors for each patch (nVect(i) = nVx(i), nVy(i), nVz(i))
    x_pc,y_pc,z_pc........ location of the point charge
    pcharge............... charge of point charge
    sigma_f............... vector of patch free charge surface densities 
    k_air................. dielectric constant of the air (or surrounding medium)
    k_obj................. dielectric constant of the sphere material

    Returns;
    sigma_b............... vector of patch bound charge surface densities
    b..................... vector of known charge (A*sigma_b = b). See Barrios and Luijten 2014, Journ. of Chem Phys. 
    

Recall: F = qE = 1/4/pi/epsilon_0*q1*q2*r/(r^3)
        E = grad(phi), phi = 1/4/pi/epsilon_0*pcharge*r/(r^2);
        Only need to consider Force of point charge on each patch. The
        patches form one body (cannot apply a force to itself)
        
        This calculation is only valid at a snapshot in time, sigma_b
        evolves over time.
%}


Npatches = length(x);
theta = acos(x./((x.^2 + y.^2).^(1/2))

%{
% Normal Vector Matrix:
% nVectM(:,:,1) = nVX1, nVX2, nVX3, ... (repeated for each row)
% nVectM(:,;,2) = nVY1, nVY2, nVY3, ... (repeated for each row)
% etc.
nVectM = zeros(Npatches, Npatches,3);

for i = 1:3
    nVectM(:,:,i) = repmat(nVect(:,i),1,Npatches);
end
%}

k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);


%{
% PATCH-PATCH INTERACTIONS NOT NEEDED WITHIN A SINGLE SPHERE! For a single sphere, just need to consider effect of external point charge on each patch

% Patch-to-patch location differences matrix
% ppld(:,:,1) = x1-x1, x1-x2, x1-x3, ... 
%               x2-x1, x2-x2, x2-x3, ... xi-xj
ppld = zeros(Npatches,Npatches,3); 
ppld(:,:,1) = x - x'; ppld(:,:,2) = y - y'; ppld(:,:,3) = z - z';
rpp = sqrt(ppld(:,:,1).*ppld(:,:,1) + ppld(:,:,2).*ppld(:,:,2) + ...
    ppld(:,:,3).*ppld(:,:,3));
%}

%%%%%%%%%%%
% NEED TO CALCULATE FORCE OF SPHERES ON EACHOTHER.
% USE SPHERE ID TO REMOVE FORCE CONTRIBUTIONS FROM PATCHES ON THE SAME SPHERE... SPHERES CAN'T APPLY FORCES TO THEMSELVES  
%%%%%%%%%%%


% pCharge-to-patch location differences vector
% pcpld(:,:) = x1-xpc, y1-ypc, z1-zpc
%              x2-xpc, y2-ypc, z2-zpc
%              etc.
pcpld = zeros(Npatches,3);
pcpld(:,1) = x-x_pc; pcpld(:,2) = y-y_pc; pcpld(:,3) = z-z_pc;
rpcp = sqrt(pcpld(:,1).*pcpld(:,1) + pcpld(:,2).*pcpld(:,2) + ...
    pcpld(:,3).*pcpld(:,3));

% nVect 
% nvx(i,:) = nvx_i, nvy_i, nvz_i (Repeat for each row i)

phi = 1/4/pi/epsilon_0 * pcharge ./rpcp;
phi_0 = 1/4/pi/epsilon_0 * pcharge ./R;
phi_norm = phi./phi_0;

end
