function [x,y,z,dA,V,nVect] = F_createSimpleCylinder(R,NpatchesCyl,dxs,dys,dzs,L)
% PROVIDES NORMALIZED POTENTIAL ENERGY OF THE POINT CHARGE SPHERE SYSTEM
%{   
    Given:
    R..................... vector of sphere radii
    NpatchesSph........... number of patches per sphere (EVENTUALLY SHOULD BE A VECTOR OF DIFFERENT VALUES FOR EACH SPHERE)
    numSpheres............ number of spheres
   
    Returns;
    x..................... vector of x-locations for each patch
    y..................... vector of y-locations for each patch
    z..................... vector of z-locations for each patch
    dA.................... vector of Area elements for each patch
    dAmat................. matrix of repeated dA vectors
    nVect................. vector of patch bound charge surface densities
    sphereID.............. vector of sphereIDs for each patch
%}


Npatches = NpatchesCyl;
x = zeros(Npatches,1)';
y = zeros(Npatches,1)';
z = zeros(Npatches,1)';
sphereID = zeros(Npatches,1);

dA = zeros(Npatches,1);
% for n = 1:numAxiSpheres
% for i = 1:NpatchesSph
%     dA(i+(n-1)*NpatchesSph) = 4*pi*(R(n)^2)/NpatchesSph;
% end
% end


% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i    
nVect = zeros(Npatches,3); % Normal Vectors

numCSs = 50; % Number of Cross Sections
numTheta = Npatches/numCSs;
theta = linspace(0,2*pi,numTheta);

for ll = 1:numCSs
for th = 1:length(theta)
    i = (ll-1)*numTheta + th;

    R_curr = R;
            
    x(i) = R_curr*cos(theta(th));
    y(i) = R_curr*sin(theta(th));
    z(i) = (ll-1)*L/numCSs;
    
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i

     nVect(i,1) = (x(i))/R_curr; 
     nVect(i,2) = (y(i))/R_curr; 
     nVect(i,3) = 0;

    % dA:
     dA(i) = R_curr*2*pi*L/Npatches; 

    % Volume:
     V = R_curr^2*pi*L;

end
end

x=x'; y=y'; z=z';

end

