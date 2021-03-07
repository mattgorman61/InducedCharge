function [x,y,z,dA,dAmat,nVect,ellID] = F_createEllipsoids(a,b,c,NpatchesEll,NEll,dxes,dyes,dzes)
%UNTITLED2 Summary of this function goes here
%{

Given:
    a..................... vector of x-radii for each ellipsoid
    b..................... vector of y-radii for each ellipsoid
    c..................... vector of z-radii for each ellipsoid
    NpatchesEll........... number of patches per sphere (EVENTUALLY SHOULD BE A VECTOR OF DIFFERENT VALUES FOR EACH SPHERE)
    NEll............ ......number of ellipsoids
   
    Returns;
    x..................... vector of x-locations for each patch
    y..................... vector of y-locations for each patch
    z..................... vector of z-locations for each patch
    dA.................... vector of Area elements for each patch
    dAmat................. matrix of repeated dA vectors
    nVect................. vector of patch bound charge surface densities
    ellID.............. vector of sphereIDs for each patch

%}

Npatches = NpatchesEll*NEll;

x = zeros(Npatches,1);
y = zeros(Npatches,1);
z = zeros(Npatches,1);
ellID = zeros(Npatches,1);

nVect = zeros(Npatches,3);


gAng = 1/2*(sqrt(5) - 1);

% Discretize Ellipsoid Surfaces
for n = 1:NEll
for i = 1:NpatchesEll
    ii = i + (n-1)*NpatchesEll;
    z(ii)= c(n)*((2*i-1)/NpatchesEll - 1) + dzes(n);
    x(ii) = a(n)*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*cos(2*pi*i*gAng) + dxes(n);
    y(ii) = b(n)*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*sin(2*pi*i*gAng) + dyes(n);
    ellID(ii) = n;
end
end

% Find nVect
for n = 1:NEll
for i = 1:NpatchesEll
    ii = i + (n-1)*NpatchesEll;
    
    zz = z(ii)-dzes(n); yy = y(ii)-dyes(n); xx = x(ii)-dxes(n);
    denom = sqrt(xx^2/a(n)^4 + yy^2/b(n)^4 + zz^2/c(n)^4);
    
    nVect(ii,3) = zz/c(n)^2/denom;
    nVect(ii,1) = xx/a(n)^2/denom;
    nVect(ii,2) = yy/b(n)^2/denom;
end
end

% Find dA
dA = [];

p = 1.6075;
A_placehold = zeros(NEll,1);

for i = 1:NEll
    A_placehold(i) = 4*pi*( ((a(i)^p).*(b(i)^p) + (a(i)^p).*(c(i)^p) + (b(i)^p).*(c(i)^p))/3 ) ^ (1/p);
end
dA_placehold = A_placehold./NpatchesEll;

for i = 1:NEll
    dA_i = dA_placehold(i).*ones(NpatchesEll,1);
    dA = [dA;dA_i];
end

% dA Matrix:
% dAmat(i) = dA(1), dA(2), ... <Npatches>
dAmat = repmat(dA',Npatches,1);




%{
figure(1);
scatter3(x,y,z,'filled');
hold on;
quiver3(x,y,z,nV(:,1),nV(:,2),nV(:,3),4,'k');
axis equal;
%}

end

