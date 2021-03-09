function [x,y,z,dA,dAmat,nVect,sphereID] = F_createSpheresEquiAng(R,NpatchesSph,numSpheres,dxs,dys,dzs)
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

Npatches = NpatchesSph*numSpheres;
x = zeros(Npatches,1)';
y = zeros(Npatches,1)';
z = zeros(Npatches,1)';
sphereID = zeros(Npatches,1);

numTheta = 40;
numPhi = NpatchesSph/numTheta;
theta = linspace(0.01,2*pi-0.01,numTheta);
phi = linspace(0.01,pi-0.01,numPhi);

%Npatches == (numTheta*numPhi*numSpheres)

dA = zeros(Npatches,1);
for n = 1:numSpheres
for i = 1:NpatchesSph
    dA(i+(n-1)*NpatchesSph) = 4*pi*(R(n)^2)/NpatchesSph;
end
end

% dA Matrix:
% dAmat(i) = dA(1), dA(2), ... <Npatches>
dAmat = repmat(dA',Npatches,1);

% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i    
nVect = zeros(Npatches,3); % Normal Vectors

%{
% Fibonacci method
gRat = (sqrt(5.0)+1.0)/2.0; % Golden Ratio
gAng = (2.0 - gRat)*(2.0*pi);
for n = 1:numSpheres
for i = 1:NpatchesSph
    lat = asin(-1.0+2.0*double(i)/(NpatchesSph+1));
    lon = gAng*i;
    
    x(i+(n-1)*NpatchesSph) = R(n)*cos(lon)*cos(lat) + dxs(n);
    y(i+(n-1)*NpatchesSph) = R(n)*sin(lon)*cos(lat) + dys(n);
    z(i+(n-1)*NpatchesSph) = R(n)*sin(lat) + dzs(n);
    sphereID(i+(n-1)*NpatchesSph) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    nVect(i+(n-1)*NpatchesSph,1) = (x(i+(n-1)*NpatchesSph)-dxs(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,2) = (y(i+(n-1)*NpatchesSph)-dys(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,3) = (z(i+(n-1)*NpatchesSph)-dzs(n))/R(n);
end
end

for n = 1:numSpheres
for i = 1:NpatchesSph
    lat = asin(-1.0+2.0*double(i)/(NpatchesSph+1));
    lon = gAng*i;
    
    x(i+(n-1)*NpatchesSph) = R(n)*cos(lon)*cos(lat) + dxs(n);
    y(i+(n-1)*NpatchesSph) = R(n)*sin(lon)*cos(lat) + dys(n);
    z(i+(n-1)*NpatchesSph) = R(n)*sin(lat) + dzs(n);
    sphereID(i+(n-1)*NpatchesSph) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    nVect(i+(n-1)*NpatchesSph,1) = (x(i+(n-1)*NpatchesSph)-dxs(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,2) = (y(i+(n-1)*NpatchesSph)-dys(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,3) = (z(i+(n-1)*NpatchesSph)-dzs(n))/R(n);
end
end

x=x'; y=y'; z=z';
%}


%% Equal Angle Method
for n = 1:numSpheres
for th = 1:length(theta)
for ph = 1:length(phi)
    
    i = (n-1)*NpatchesSph + (th-1)*length(phi) + ph;
    
    xx = R*sin(phi(ph))*cos(theta(th));% + dxs(n);
    yy = R*sin(phi(ph))*sin(theta(th));% + dys(n);
    zz = R*cos(phi(ph));% + dzs(n);
    
    x(i) = xx;
    y(i) = yy;
    z(i) = zz;
    %nVect(i,1) = (x(i)-dxs(n))/R(n); 
    %nVect(i,2) = (y(i)-dys(n))/R(n); 
    %nVect(i,3) = (z(i)-dzs(n))/R(n);
    
    nVect(i,1) = (xx/R); 
    nVect(i,2) = (yy/R); 
    nVect(i,3) = (zz/R);
    
end
end
end
%}

%{
% Length Check:
length(nVect(:,3))
length(z)
%}

x=x'; y=y'; z=z';
    

end

