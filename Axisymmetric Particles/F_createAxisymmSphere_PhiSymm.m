function [x,y,z,dA,nVect,sphereID] = F_createAxisymmSphere_PhiSymm(R,NpatchesSph,dxs,dys,dzs,n)
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

lPhiSymm = true; % Logical: determines whether symmetry is about phi or about theta

Npatches = NpatchesSph;
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

% Axisymmetric Particles
numTheta = 50;

theta = linspace(0,2*pi,numTheta);
numPhi = NpatchesSph/numTheta;
phi = linspace(0,pi,numPhi);

for th = 1:numTheta
for ph = 1:numPhi
    i = (th-1)*numPhi + ph;

    R_curr = R + R*0.1*sin(10*phi(ph)); % Opposite trig operations for stronger attraction, same for repulsion?
%   R_curr = 1;
            
    x(i) = R_curr*cos(theta(th))*sin(phi(ph)) + dxs;
    y(i) = R_curr*sin(theta(th))*sin(phi(ph)) + dys;
    z(i) = R_curr*cos(phi(ph)) + dzs;
    sphereID(i) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
     nVect(i,1) = (x(i)-dxs)/R_curr; 
     nVect(i,2) = (y(i)-dys)/R_curr; 
     nVect(i,3) = (z(i)-dzs)/R_curr;
     
    % dA:
     dA(i) = R_curr^2*2*pi/numTheta*pi/numPhi; 

end
end

x=x'; y=y'; z=z';

end

