function [x,y,z,dA,dAmat,nVect,ellID,a,b,c,x0,y0,z0,psi,phi,theta,A_Eul,EllPatchData] = F_createEllipsoid(a,b,c,NpatchesEll,NEll,EllPatchData,dx,dy,dz,psi,phi,theta,ellID_in)
% Creates Geometry for a single ellipsoid
%{

Given:
    a..................... vector of x-radii for each ellipsoid
    b..................... vector of y-radii for each ellipsoid
    c..................... vector of z-radii for each ellipsoid
    NpatchesEll........... number of patches per sphere (EVENTUALLY SHOULD BE A VECTOR OF DIFFERENT VALUES FOR EACH SPHERE)
    NEll.................. number of ellipsoids
    EllpatchData.......... matrix with Ellipsoid patch data
    dx.................... x location of the ellipsoid center
    dy.................... y location of the ellipsoid center
    dz.................... z location of the ellipsoid center
    psi................... Euler rotation angle
    phi................... Euler rotation angle
    theta................. Euler rotation angle
    EllID................. EllipseID

   
Returns:
    x..................... vector of x-locations for each patch
    y..................... vector of y-locations for each patch
    z..................... vector of z-locations for each patch
    dA.................... vector of Area elements for each patch
    dAmat................. matrix of repeated dA vectors
    nVect................. vector of patch bound charge surface densities
    ellID................. vector of particle IDs for each patch

Updates:
    EllpatchData

%}

Npatches = NpatchesEll*NEll;

x1 = zeros(NpatchesEll,1);
y1 = zeros(NpatchesEll,1);
z1 = zeros(NpatchesEll,1);
nVect = zeros(NpatchesEll,3);

x = zeros(NpatchesEll,1);
y = zeros(NpatchesEll,1);
z = zeros(NpatchesEll,1);
ellID = zeros(NpatchesEll,1);



% Discretize Ellipsoid Surfaces, centered at origin
gAng = 1/2*(sqrt(5) - 1);
for i = 1:NpatchesEll
    % ii = i + (n-1)*NpatchesEll;
    z1(i)= c*((2*i-1)/NpatchesEll - 1);
    x1(i) = a*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*cos(2*pi*i*gAng);
    y1(i) = b*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*sin(2*pi*i*gAng);
    ellID(i) = ellID_in;
end


% Apply rotation
for n = 1:NEll
    [xx,yy,zz,A_Eul_dim3] = F_EulRot(x1,y1,z1,psi,phi,theta ); % Rotate the Ellipse
    x = xx; y = yy; z = zz;
    A_Eul = [A_Eul_dim3, [0 0 0]']; % Have to use 4x4 Rotation Matrix, where last row = [0 0 0 1].
    A_Eul = [A_Eul; [0 0 0 1]];
end


% Apply offset
for i = 1:NpatchesEll
    % ii = i + (n-1)*NpatchesEll;
    z(i)= z(i) + dz;
    x(i) = x(i) + dx;
    y(i) = y(i) + dy;
end


% Find nVect
for i = 1:NpatchesEll    
    zz = z(i)-dz; yy = y(i)-dy; xx = x(i)-dx;
    denom = sqrt(xx^2/a^4 + yy^2/b^4 + zz^2/c^4);
    
    nVect(i,3) = zz/c^2/denom;
    nVect(i,1) = xx/a^2/denom;
    nVect(i,2) = yy/b^2/denom;
end


% Find dA
p = 1.6075;
A_ellipse = 4*pi*( ((a^p).*(b^p) + (a^p).*(c^p) + (b^p).*(c^p))/3 ) ^ (1/p); % Area of the Ellipse
dA_ellipse = A_ellipse/NpatchesEll; % Area of each patch

dA = ones(NpatchesEll,1).*dA_ellipse; % Vector of patch area

% dA Matrix:
% dAmat(i) = dA(1), dA(2), ... <Npatches>
dAmat = repmat(dA',Npatches,1);

x0 = dx;
y0 = dy;
z0 = dz;   


n = ellID_in;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),1) = x;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),2) = y;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),3) = z;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),4) = dA;

EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),5) = nVect(:,1);
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),6) = nVect(:,2);
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),7) = nVect(:,3);

EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),8) = a;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),9) = b;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),10) = c;

EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),11) = x0;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),12) = y0;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),13) = z0;    

EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),14) = psi;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),15) = phi;
EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),16) = theta;

EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),17) = ellID;





end

