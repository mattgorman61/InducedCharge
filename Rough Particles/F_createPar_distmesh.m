function [x,y,z,dA,nVect,sphereID] = F_createEllipPar_distmesh(R,NpatchesSph,dxs,dys,dzs,n)

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

a = 2; b = 1; c = 1.5;
R_0 = 1;
alph = 5;
A = 0.15;

% [xx,yy,zz] = meshgrid(-R_0:0.1:R_0, -R_0:0.1:R_0, -R_0:0.1:R_0);
% dd = sqrt(xx.^2 + yy.^2 + zz.^2) - R_0.*(1 + A.*sin(alph*atan2(sqrt(xx.^2 + yy.^2),zz)));
% [p,t]=distmeshsurface(@dmatrix,@huniform,0.2,[-R_0,-R_0,-R_0; R_0,R_0,R_0],[],xx,yy,zz);

 fd=@(p) sqrt(p(:,1).^2 + p(:,2).^2 + p(:,3).^2) - R_0*(1+A*sin(alph*atan2(sqrt(p(:,1).^2 + p(:,2).^2),p(:,3)))  );
 [p,t]=distmeshsurface(fd,@huniform,0.25,1.5*[-R_0,-R_0,-R_0; R_0,R_0,R_0]);

x = p(:,1); 
y = p(:,2);
z = p(:,3);

figure;
axis equal; box on; 
scatter3(x,y,z,'filled');
axis equal;


dA = 0;
nVect = 0;
sphereID = n;


end






















































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
