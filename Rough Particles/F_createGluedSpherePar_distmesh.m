function [x,y,z,dA,V,nVect] = F_createGluedSpherePar_distmesh(R,dxs,dys,dzs,elemSize)
% PROVIDES ...
%{  

  [P,T]=DISTMESHSURFACE(FD,FH,H0,BBOX,FPARAMS)

  NOTE:
     Doesn't work too well for nonuniform size functions... need a
     better initial grid or density control. Also no support for
     fixed points, again because of the initial grid.

     P:         Node positions (Nx2)
     T:         Triangle indices (NTx3)
     FD:        Distance function d(x,y)
     FH:        Scaled edge length function h(x,y)
     H0:        Initial edge length
     BBOX:      Bounding box [xmin,ymin; xmax,ymax]
     FPARAMS:   Additional parameters passed to FD and FH

Example: (Uniform Mesh on Ellipsoid)
     fd=@(p) p(:,1).^2/4+p(:,2).^2/1+p(:,3).^2/1.5^2-1;
     [p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);

%}

currFolder = pwd;
% fprintf('%s',currFolder);
path1 = strcat(currFolder,'\..\Functions');        addpath(path1);
path2 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\distmesh\distmesh'; addpath(path2);

% First step: save sphere center locations to current directory
sphData = [dxs',dys',dzs',R'];
dirString = pwd;
fileName = [dirString,'\sph_Data.mat'];
save(fileName,'sphData')

R_0 = R;


[p,t]=distmeshsurface_gluedPars(@fd_gluedSpheres,@huniform,elemSize,7*[-R_0,-R_0,-R_0; R_0,R_0,R_0]);

x = p(:,1); 
y = p(:,2);
z = p(:,3);

r_spheres = sqrt(dxs.^2 + dys.^2 + dzs.^2);

figure;
axis equal; box on; 
scatter3(x,y,z,'filled');
axis equal;

NpatchesSph = length(x);
dA = zeros(NpatchesSph,1);
nVect = zeros(NpatchesSph,3);

for i = 1:NpatchesSph
    
    x_i = x(i); y_i = y(i); z_i = z(i);
    r_candidates = zeros(length(dxs),1);
    for n = 1:length(r_candidates)
        r_candidates(n) = sqrt((x_i-dxs(n)).^2 + (y_i-dys(n)).^2 + (z_i-dzs(n)).^2 );
    end

    A = abs(r_candidates - R_0);
    minRadDiff = min(min(A));
    [sphInd,rCandInd] = find(A==minRadDiff);
    minRadDiffCheck = A(rCandInd,sphInd);
    
%     r_i = r_candidates(rCandInd);
    sphereXLoc = dxs(sphInd(1));
    sphereYLoc = dys(sphInd(1));
    sphereZLoc = dzs(sphInd(1));
    r_i = sqrt((x_i-sphereXLoc)^2 + (y_i-sphereYLoc)^2 + (z_i-sphereZLoc)^2);

     nVect(i,1) = (x_i-sphereXLoc)/r_i; 
     nVect(i,2) = (y_i-sphereYLoc)/r_i; 
     nVect(i,3) = (z_i-sphereZLoc)/r_i;

end

shp = alphaShape(x,y,z);
figure;
plot(shp);

totalSurfArea = surfaceArea(shp);
dA = totalSurfArea/length(x);
V = volume(shp);


end

function d=fd_gluedSpheres(p)

    
    load('sph_Data.mat');
    sphLocX = sphData(:,1);
    sphLocY = sphData(:,2);
    sphLocZ = sphData(:,3);
    sphRad = sphData(:,4);
    
    d = sqrt( (p(:,1)-sphLocX(1)).^2 + (p(:,2)-sphLocY(1)).^2 + (p(:,3)-sphLocZ(1)).^2) - sphRad(1)^2;
    
    for i = 2:length(sphLocX)
        d_i = sqrt( (p(:,1)-sphLocX(i)).^2 + (p(:,2)-sphLocY(i)).^2 + (p(:,3)-sphLocZ(i)).^2) - sphRad(i)^2;
        d = dunion(d,d_i);
    end
    
end









    % dA
%     A3 =  A^2*(1 + (-cos(2*alpha*pi)-2*alpha)/2/(1+2*alpha));
%     A2 = 2*A*sin(alpha*pi);
%     A1 = 2;
%     SA = 2*pi*(R_0^2)*(A1 + A2 + A3);
%     
%     dA = SA/NpatchesSph*ones(length(x),1); 
%     
%     
%     V3 = A^3/(1+9*alpha^2)*((sin(alpha*pi)^3) - 18*alpha^2);
%     V2 = 3*A^2/(1-4*alpha^2)*((sin(alpha*pi)^2) - 4*alpha^2);
%     V1 = 3*A/(1-alpha^2)*(sin(alpha*pi));
%     V = 2*pi*(R_0^3)/3*(2 + V1 + V2 + V3);











































%{
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
buffer = 0.15;
RippleAmp = 0.2;
alpha = 4;
numTheta = NpatchesSph/20;
theta = linspace(0+buffer,2*pi-buffer,numTheta);
bias = 1;



numPhi = NpatchesSph/numTheta;
phi = [F_createVector_bias(1/bias,numPhi/2,pi/2 - buffer) + buffer; pi/2 + F_createVector_bias(bias,numPhi/2,pi/2 - buffer)];
% phi = [F_createVector_bias(1/bias,numPhi/2,pi/2); pi/2 + F_createVector_bias(bias,numPhi/2,pi/2)];

figure;
scatter(phi,zeros(length(phi),1));

for th = 1:numTheta
for ph = 1:numPhi
    i = (th-1)*numPhi + ph;

    R_curr = R + R*RippleAmp*sin(alpha*phi(ph)); 
%   R_curr = 1;
            
    x(i) = R_curr*cos(theta(th))*sin(phi(ph)) + dxs;
    y(i) = R_curr*sin(theta(th))*sin(phi(ph)) + dys;
    z(i) = R_curr*cos(phi(ph)) + dzs;
    sphereID(i) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    dPhi_r = R*(RippleAmp*alpha*cos(alpha*(phi(ph))));
    
    nvx = R_curr*sin(phi(ph))*cos(theta(th))*(dPhi_r*cos(phi(ph)) - R_curr*sin(phi(ph)));
    nvy = R_curr*sin(phi(ph))*sin(theta(th))*(dPhi_r*cos(phi(ph)) - R_curr*sin(phi(ph)));
    nvz = -R_curr*sin(phi(ph))*sin(theta(th))*(dPhi_r*sin(phi(ph)) + R_curr*cos(phi(ph)))*sin(theta(th)) + ...
            - R_curr*sin(phi(ph))*cos(theta(th))*(dPhi_r*sin(phi(ph)) + R_curr*cos(phi(ph)))*cos(theta(th));
        
    nvMag = sqrt(nvx*nvx + nvy*nvy + nvz*nvz);
    if(nvMag == 0)
        nvMag = 1;
    end
    
    
    % INCORRECT!!!!! ONLY FOR SPHERES
%      nVect(i,1) = (x(i)-dxs)/R_curr; 
%      nVect(i,2) = (y(i)-dys)/R_curr; 
%      nVect(i,3) = (z(i)-dzs)/R_curr;

    if(nvMag == 0)
        nVect(i,1) = (x(i)-dxs)/R_curr; 
        nVect(i,2) = (y(i)-dys)/R_curr; 
        nVect(i,3) = (z(i)-dzs)/R_curr;
    else
        nVect(i,1) = -nvx/real(nvMag);
        nVect(i,2) = -nvy/real(nvMag);
        nVect(i,3) = -nvz/real(nvMag);    
    end

    % dA:
     dA(i) = R_curr*2*pi/numTheta * R_curr*pi/numPhi; 

end
end

x=x'; y=y'; z=z';

end
%}
