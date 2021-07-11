function [x,y,z,dA,V,nVect] = F_createAxiSymmPar_distmesh(R,dxs,dys,dzs)
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

R_0 = 1;
alpha = 15;
A = 0.15;
elemSize = 0.12;

% [xx,yy,zz] = meshgrid(-R_0:0.1:R_0, -R_0:0.1:R_0, -R_0:0.1:R_0);
% dd = sqrt(xx.^2 + yy.^2 + zz.^2) - R_0.*(1 + A.*sin(alph*atan2(sqrt(xx.^2 + yy.^2),zz)));
% [p,t]=distmeshsurface(@dmatrix,@huniform,0.2,[-R_0,-R_0,-R_0; R_0,R_0,R_0],[],xx,yy,zz);


fd=@(p) sqrt(p(:,1).^2 + p(:,2).^2 + p(:,3).^2) - R_0*(1+A*sin(alpha*atan2(sqrt(p(:,1).^2 + p(:,2).^2),p(:,3)))  );
[p,t]=distmeshsurface(fd,@huniform,elemSize,1.5*[-R_0,-R_0,-R_0; R_0,R_0,R_0]);

x = p(:,1); 
y = p(:,2);
z = p(:,3);

figure;
axis equal; box on; 
scatter3(x,y,z,'filled');
axis equal;

NpatchesSph = length(x);
dA = zeros(NpatchesSph,1);
nVect = zeros(NpatchesSph,1);

for i = 1:NpatchesSph
    
    theta_i = atan2(y(i),x(i));
    phi_i = atan2(sqrt(x(i)^2 + y(i)^2),z(i) );
    R_curr = R_0 + R_0*A*sin(alpha*phi_i); 
    %   R_curr = 1;
            
    x(i) = R_curr*cos(theta_i)*sin(phi_i) + dxs;
    y(i) = R_curr*sin(theta_i)*sin(phi_i) + dys;
    z(i) = R_curr*cos(phi_i) + dzs;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    dPhi_r = R_0*(A*alpha*cos(alpha*(phi_i)));
    
    nvx = R_curr*sin(phi_i)*cos(theta_i)*(dPhi_r*cos(phi_i) - R_curr*sin(phi_i));
    nvy = R_curr*sin(phi_i)*sin(theta_i)*(dPhi_r*cos(phi_i) - R_curr*sin(phi_i));
    nvz = -R_curr*sin(phi_i)*sin(theta_i)*(dPhi_r*sin(phi_i) + R_curr*cos(phi_i))*sin(theta_i) + ...
            - R_curr*sin(phi_i)*cos(theta_i)*(dPhi_r*sin(phi_i) + R_curr*cos(phi_i))*cos(theta_i);
        
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

end

    % dA
    A3 =  A^2*(1 + (-cos(2*alpha*pi)-2*alpha)/2/(1+2*alpha));
    A2 = 2*A*sin(alpha*pi);
    A1 = 2;
    SA = 2*pi*(R_0^2)*(A1 + A2 + A3);
    
    dA = SA/NpatchesSph*ones(length(x),1); 
    
    
    V3 = A^3/(1+9*alpha^2)*((sin(alpha*pi)^3) - 18*alpha^2);
    V2 = 3*A^2/(1-4*alpha^2)*((sin(alpha*pi)^2) - 4*alpha^2);
    V1 = 3*A/(1-alpha^2)*(sin(alpha*pi));
    V = 2*pi*(R_0^3)/3*(2 + V1 + V2 + V3);
    
    
    

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
