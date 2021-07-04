function [Fnet,Fx,Fy,Fz,F0] = F_getForces_Mult_Matrix(numPars,NpatchesPar,R,x,y,z,dA,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0)
% PROVIDES VECTOR OF FORCES ACTING ON EACH PATCH OF THE SPHERES
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
    

Recall: F = qE = 1/4/pi/epsilon_0*q1*q2*r/(r^3)
        Only need to consider Force of point charge on each patch. The
        patches form one body (cannot apply a force to itself)
        
        This calculation is only valid at a snapshot in time, sigma_b
        evolves over time.
%}

Npatches = length(x);
% dA = 4*pi*(R^2)/Npatches;


% Normal Vector Matrix:
% nVectM(:,:,1) = nVX1, nVX2, nVX3, ... (repeated for each row)
% nVectM(:,;,2) = nVY1, nVY2, nVY3, ... (repeated for each row)
% etc.
nVectM = zeros(Npatches, Npatches,3);
for i = 1:3
    nVectM(:,:,i) = repmat(nVect(:,i),1,Npatches);
end

k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);


%%{
% NOT NEEDED FOR A SINGLE PARTICLE! Just need to consider effect of point charge on each patch
% Patch-to-patch location differences matrix
% ppld(:,:,1) = x1-x1, x1-x2, x1-x3, ... 
%               x2-x1, x2-x2, x2-x3, ... xi-xj
ppld = zeros(Npatches,Npatches,3); 
ppld(:,:,1) = x - x'; ppld(:,:,2) = y - y'; ppld(:,:,3) = z - z';
rpp = sqrt(ppld(:,:,1).*ppld(:,:,1) + ppld(:,:,2).*ppld(:,:,2) + ...
    ppld(:,:,3).*ppld(:,:,3));
%}

% CHECK HERE
rpp(rpp<1) = 1; % Don't want to divide 0 by 0. These terms should go to 0 anyways due to 0 in numerator.

% Neglect contributions from patches belonging to the same particle
for i = 1:numPars
    ppld((i-1)*NpatchesPar+1:i*NpatchesPar,(i-1)*NpatchesPar+1:i*NpatchesPar,1) = 0;
    ppld((i-1)*NpatchesPar+1:i*NpatchesPar,(i-1)*NpatchesPar+1:i*NpatchesPar,2) = 0;
    ppld((i-1)*NpatchesPar+1:i*NpatchesPar,(i-1)*NpatchesPar+1:i*NpatchesPar,3) = 0;
%     rpp((i-1)*NpatchesSph+1:i*NpatchesSph,(i-1)*NpatchesSph+1:i*NpatchesSph) = 0;
end

chargevect = sigma.*dA;

chgvectSq = chargevect*chargevect';

fxTemp = ppld(:,1)./(rpp.^3).*chgvectSq * 1/4/pi/epsilon_0;
fxTemp(isnan(fxTemp) | isinf(fxTemp) | rpp < 1) = 0;
% fxTemp(isnan(fxTemp) | isinf(fxTemp) ) = 0;
Fx = sum(fxTemp)';

fyTemp = ppld(:,2)./(rpp.^3).*chgvectSq * 1/4/pi/epsilon_0;
fyTemp(isnan(fyTemp) | isinf(fyTemp) | rpp < 1) = 0;
% fyTemp(isnan(fxTemp) | isinf(fxTemp) ) = 0;
Fy = sum(fyTemp)';

fzTemp = ppld(:,3)./(rpp.^3).*chgvectSq * 1/4/pi/epsilon_0;
fzTemp(isnan(fzTemp) | isinf(fzTemp) | rpp < 1) = 0;
% fzTemp(isnan(fxTemp) | isinf(fxTemp) ) = 0;
Fz = sum(fzTemp)';

% pCharge-to-patch location differences vector
% pcpld(:,:) = x1-xpc, y1-ypc, z1-zpc
%              x2-xpc, y2-ypc, z2-zpc
%              etc.
pcpld = zeros(Npatches,3);
pcpld(:,1) = x-x_pc; pcpld(:,2) = y-y_pc; pcpld(:,3) = z-z_pc;
rpcp = sqrt(pcpld(:,1).*pcpld(:,1) + pcpld(:,2).*pcpld(:,2) + ...
    pcpld(:,3).*pcpld(:,3));


Fx = Fx + pcpld(:,1)./(rpcp.^3) * 1/4/pi/epsilon_0 * pcharge .*sigma.*dA;
Fy = Fy + pcpld(:,2)./(rpcp.^3) * 1/4/pi/epsilon_0 * pcharge .*sigma.*dA;
Fz = Fz + pcpld(:,3)./(rpcp.^3) * 1/4/pi/epsilon_0 * pcharge .*sigma.*dA;

Fnet = zeros(numPars,3);
for i = 1:numPars
    Fnet(i,:) = [sum(Fx((i-1)*NpatchesPar+1:i*NpatchesPar)),sum(Fy((i-1)*NpatchesPar+1:i*NpatchesPar)),sum(Fz((i-1)*NpatchesPar+1:i*NpatchesPar))];
end

if (pcharge == 0)
    F0 = (100^2)/16/pi/epsilon_0/(R^2);
else
    F0 = (pcharge^2)/16/pi/epsilon_0/(R^2);
end

end
