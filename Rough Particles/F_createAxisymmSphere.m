function [x,y,z,dA,Volume,nVect,sphereID] = F_createAxisymmSphere(R,NpatchesSph,dxs,dys,dzs,n)
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
RippleAmpPh = 0.1;
RippleAmpTh = 0.1;
alpha = 5;
beta = 15;
numTheta = NpatchesSph/100;
theta = linspace(0,2*pi,numTheta);
bias = 5;

numPhi = NpatchesSph/numTheta;
phi = [F_createVector_bias(1/bias,numPhi/2,pi/2); pi/2 + F_createVector_bias(bias,numPhi/2,pi/2)];
% phi = [F_createVector_bias(1/bias,numPhi/2,pi/2); pi/2 + F_createVector_bias(bias,numPhi/2,pi/2)];

% figure;
% scatter(phi,zeros(length(phi),1));

for th = 1:numTheta
for ph = 1:numPhi
    i = (th-1)*numPhi + ph;

    R_curr = R*(1 + RippleAmpPh*sin(alpha*phi(ph)))*(1 + RippleAmpTh*sin(beta*theta(th))); 
%   R_curr = 1;
            
    x(i) = R_curr*cos(theta(th))*sin(phi(ph)) + dxs;
    y(i) = R_curr*sin(theta(th))*sin(phi(ph)) + dys;
    z(i) = R_curr*cos(phi(ph)) + dzs;
    sphereID(i) = n;
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    dPhi_r = R*(RippleAmpPh*alpha*cos(alpha*(phi(ph))));
    
%     nvx = R_curr*sin(phi(ph))*cos(theta(th))*(dPhi_r*cos(phi(ph)) - R_curr*sin(phi(ph)));
%     nvy = R_curr*sin(phi(ph))*sin(theta(th))*(dPhi_r*cos(phi(ph)) - R_curr*sin(phi(ph)));
%     nvz = -R_curr*sin(phi(ph))*sin(theta(th))*(dPhi_r*sin(phi(ph)) + R_curr*cos(phi(ph)))*sin(theta(th)) + ...
%             - R_curr*sin(phi(ph))*cos(theta(th))*(dPhi_r*sin(phi(ph)) + R_curr*cos(phi(ph)))*cos(theta(th));
%         
%     nvMag = sqrt(nvx*nvx + nvy*nvy + nvz*nvz);
%     if(nvMag == 0)
%         nvMag = 1;
%     end
    
    
    % INCORRECT!!!!! ONLY FOR SPHERES
     nVect(i,1) = (x(i)-dxs)/R_curr; 
     nVect(i,2) = (y(i)-dys)/R_curr; 
     nVect(i,3) = (z(i)-dzs)/R_curr;

%     nVect(i,1) = -nvx/real(nvMag);
%     nVect(i,2) = -nvy/real(nvMag);
%     nVect(i,3) = -nvz/real(nvMag);     
     
    % dA:
     dA(i) = R_curr*2*pi/numTheta * R_curr*pi/numPhi; 

end
end

x=x'; y=y'; z=z';
Volume = NaN;

end

