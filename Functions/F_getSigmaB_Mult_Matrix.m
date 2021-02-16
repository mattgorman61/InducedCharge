function [sigma_b2,b2] = F_getSigmaB_Mult_Matrix(numSpheres,NpatchesSph,R,dA,dAmat,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z)
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
    b..................... vector of known charge (A*sigma_b = b). See Barros and Luijten 2014, Journ. of Chem Phys. 
    
%}

Npatches = length(x);

% Normal Vector Matrix:
% nVectM(:,:,1) = nVX1, nVX2, nVX3, ... (repeated for each row)
% nVectM(:,;,2) = nVY1, nVY2, nVY3, ... (repeated for each row)
% etc.
nVectM = zeros(Npatches, Npatches,3);
for i = 1:3
nVectM(:,:,i) = repmat(nVect(:,i),1,Npatches);
end

k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);

% Patch-to-patch location differences matrix
% ppld(:,:,1) = x1-x1, x1-x2, x1-x3, ... 
%               x2-x1, x2-x2, x2-x3, ... xi-xj
ppld = zeros(Npatches,Npatches,3); 
ppld(:,:,1) = x - x'; ppld(:,:,2) = y - y'; ppld(:,:,3) = z - z';
rpp = sqrt(ppld(:,:,1).*ppld(:,:,1) + ppld(:,:,2).*ppld(:,:,2) + ...
    ppld(:,:,3).*ppld(:,:,3));

% Dot product term of A
M = (ppld(:,:,1).*nVectM(:,:,1) + ppld(:,:,2).*nVectM(:,:,2) + ...
    ppld(:,:,3).*nVectM(:,:,3))./(rpp.^3);
M(isnan(M) | isinf(M)) = 0;

A2 = k_bar*eye(Npatches) + k_delta/4/pi*dAmat.*M;

% pCharge-to-patch location differences vector
% pcpld(:,:) = x1-xpc, y1-ypc, z1-zpc
%              x2-xpc, y2-ypc, z2-zpc
%              etc.
pcpld = zeros(Npatches,3);
pcpld(:,1) = x-x_pcs; pcpld(:,2) = y-y_pcs; pcpld(:,3) = z-z_pcs;
rpcp = sqrt(pcpld(:,1).*pcpld(:,1) + pcpld(:,2).*pcpld(:,2) + ...
    pcpld(:,3).*pcpld(:,3));

% nVect 
% nvx(i,:) = nvx_i, nvy_i, nvz_i (Repeat for each row i)

% Pcharge-to-patch dot product matrix
Mpc = diag(pcpld*nVect')./(rpcp.^3);
%Mpc = 
Mpc(isnan(Mpc) | isinf(Mpc)) = 0;

b2 = ((1-k_bar)*eye(Npatches) - k_delta/4/pi*dAmat.*M)*sigma_f - k_delta/4/pi*pcharge*Mpc - ...
    k_delta/4/pi*(Ext_EField_x.*nVect(:,1)+ Ext_EField_y.*nVect(:,2) + Ext_EField_z.*nVect(:,3));


sigma_b2 = A2\b2;
