function [x,y,z,dA,dAmat,nVect,sphereID] = F_createSingSphere(R,NpatchesSph,dxs,dys,dzs,n)
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

x = zeros(NpatchesSph,1)';
y = zeros(NpatchesSph,1)';
z = zeros(NpatchesSph,1)';
sphereID = zeros(NpatchesSph,1);

dA = zeros(NpatchesSph,1);

for i = 1:NpatchesSph
    dA(i) = 4*pi*(R^2)/NpatchesSph;
end


% dA Matrix:
% dAmat(i) = dA(1), dA(2), ... <Npatches>
dAmat = repmat(dA',NpatchesSph,1);

% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i    
nVect = zeros(NpatchesSph,3); % Normal Vectors

% Fibonacci method
gRat = (sqrt(5.0)+1.0)/2.0; % Golden Ratio
gAng = (2.0 - gRat)*(2.0*pi);



for i = 1:NpatchesSph
    lat = asin(-1.0+2.0*double(i)/(NpatchesSph+1));
    lon = gAng*i;
    
    x(i) = R*cos(lon)*cos(lat) + dxs;
    y(i) = R*sin(lon)*cos(lat) + dys;
    z(i) = R*sin(lat) + dzs;
    sphereID(i) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    nVect(i,1) = (x(i)-dxs)/R; 
    nVect(i,2) = (y(i)-dys)/R; 
    nVect(i,3) = (z(i)-dzs)/R;
    

end

x=x'; y=y'; z=z';

end

