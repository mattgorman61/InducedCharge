clc; clear all; close all;

%% Create Multiple Ellipsoids

a = [1 1 1]; b = [2 2 2]; c = [1 1 1];
dxes = [1 3 5]; dyes = [1 3 4]; dzes = [0 0 0];
psi = [0 pi/4 pi/2]; % Rotation about z-axis
phi = [0 0 0]; % Rotation about x-axis
theta = [0 0 0]; % Rotation about y-axis
NpatchesEll = 1000; 
NEll = 3;
Npatches = NEll*NpatchesEll;
EllPatchData = zeros(Npatches,14); % X,Y,Z,dA, nVx,nVy,nVz, a,b,c, x0,y0,z0, ellID of each patch

tsteps = 100;

for n = 1:NEll
    [x,y,z,dA,dAmat,nVect,ellID,a1,b1,c1,x0,y0,z0,EllPatchData] = F_createEllipsoid(a(n),b(n),c(n),NpatchesEll,NEll,EllPatchData,dxes(n),dyes(n),dzes(n),psi(n),phi(n),theta(n),n);
    
    %{
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),1) = x;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),2) = y;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),3) = z;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),4) = dA;
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),5) = nVect(:,1);
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),6) = nVect(:,2);
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),7) = nVect(:,3);
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),8) = a1;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),9) = b1;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),10) = c1;
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),11) = x0;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),12) = y0;
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),13) = z0;    
    
    EllPatchData((n + (n-1)*NpatchesEll):(n*(1+NpatchesEll)-1),14) = ellID;
    %}
end

colors = ['r','b','k','g','m'];
figure();
hold on;
for n = 1:NEll
    x_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),1);
    y_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),2); 
    z_i = EllPatchData((1+(n-1)*NpatchesEll):((n)*NpatchesEll),3);
    scatter3(x_i,y_i,z_i,'filled',colors(n));
    % quiver3(x_i,y_i,z_i,nVect(:,1),nVect(:,2),nVect(:,3),4); % Show Normal Vectors
end
axis equal;
grid on;
view(45,25);
xlabel('x'); ylabel('y'); zlabel('z');

%for i=1:tsteps
collisionList = F_CollideCheckCall(EllPatchData, NEll, NpatchesEll)
%end

%% Test Section

%%{
% BEGIN TEST SECTION
fprintf('\n\n BEGIN TEST SECTION: \n\n');

a0 = 1; b0 = 2; c0 = 1;
x0 = 1; y0 = 1; z0 = 1;
a1 = 1; b1 = 2; c1 = 1;
x1 = 1; y1 = 2; z1 = 1;

lCollideFlag = -1;

A0 = 1/a0^2; B0 = 1/b0^2; C0 = 1/c0^2; 
D0 = 0; E0 = 0; F0 = 0;
G0 = -2*x0/a0^2; H0 = -2*y0/b0^2; J0 = -2*z0/c0^2;
K0 = (x0^2)/a0^2 + (y0^2)/b0^2 + (z0^2)/c0^2 - 1;

S_0 = [2*A0, D0, F0, G0; D0, 2*B0, E0, H0; F0, E0, 2*C0, J0; G0,H0,J0,2*K0];
    
A1 = 1/a1^2; B1 = 1/b1^2; C1 = 1/c1^2; 
D1 = 0; E1 = 0; F1 = 0;
G1 = -2*x1/a1^2; H1 = -2*y1/b1^2; J1 = -2*z1/c1^2;
K1 = (x1^2)/a1^2 + (y1^2)/b1^2 + (z1^2)/c1^2 - 1;

S_1 = [2*A1, D1, F1, G1; D1, 2*B1, E1, H1; F1, E1, 2*C1, J1; G1,H1,J1,2*K1];

EigMat = (S_0)\S_1;
[EigVect,~] = eig(EigMat);

lCollideFlag = -1;
for i=1:length(EigVect(1,:))
    evecti = EigVect(:,i);
    if(isequal(evecti,[x0;y0;z0;1]))
        lCollideFlag = 1;
    end
    if(~isequal(imag(evecti),[0;0;0;0]))
        lCollideFlag = 1;
    end
end

lCollideFlag

% END TEST SECTION
%} 

%% Plot Original Ellipse
%{
NpatchesEll = 200;
x = zeros(NpatchesEll,1); y = zeros(NpatchesEll,1); z = zeros(NpatchesEll,1);

n = 1;
a = 4; b = 1; c = 1;
gAng = 1/2*(sqrt(5) - 1);
for i = 1:NpatchesEll
    ii = i + (n-1)*NpatchesEll;
    z(ii) = c(n)*((2*i-1)/NpatchesEll - 1);
    x(ii) = a(n)*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*cos(2*pi*i*gAng);
    y(ii) = b(n)*sqrt( (2*i-1)/NpatchesEll*(2-(2*i-1)/NpatchesEll))*sin(2*pi*i*gAng);
    %ellID(ii) = n;
end

figure();
scatter3(x,y,z);

[xrot,yrot,zrot] = F_EulRot(x,y,z,pi/4,pi/4,pi/4);
%[xrot,yrot,zrot] = F_QuatRot(x,y,z,pi/4,pi/4,pi/4);
hold on;
scatter3(xrot,yrot,zrot,'filled','r');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');

%}


%% Euler Angle Matrix: Axis Rotation

%{
psi = -pi/4;
phi = -pi/4;
theta = pi/4;

A_Eul = zeros(3);

A_Eul(1,1) = cos(psi)*cos(phi) - cos(theta)*sin(theta)*sin(psi);
A_Eul(1,2) = -sin(psi)*cos(theta) - cos(theta)*sin(phi)*cos(psi);
A_Eul(1,3) = sin(theta)*sin(phi);

A_Eul(2,1) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
A_Eul(2,2) = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
A_Eul(2,3) = -sin(theta)*cos(phi);

A_Eul(3,1) = sin(theta)*sin(psi);
A_Eul(3,2) = sin(theta)*cos(psi);
A_Eul(3,3) = cos(theta);

fig = figure();
x = [0 0 0]; y = [0 0 0]; z = [0 0 0];
%x = x'; y = y'; z = z';
vx = [1 0 0]; vy = [0 1 0]; vz = [0 0 1];
%vx = vx'; vy = vy'; vz = vz;
quiver3(x,y,z,vx,vy,vz,'k');

vx_rot = A_Eul*vx'; vy_rot = A_Eul*vy'; vz_rot = A_Eul*vz';
hold on;
quiver3(x,y,z,vx_rot',vy_rot',vz_rot','r');
axis equal;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
xlabel('x'); ylabel('y'); zlabel('z');

%% GET QUATERNIONS

%phi = -phi;
%psi = -psi;
%theta = -theta;

%Definition from Adhesive Particle Flow
%%{
eps1 = cos((phi-psi)/2)*sin(theta/2); 
eps2 = sin((phi-psi)/2)*sin(theta/2);
eps3 = sin((phi+psi)/2)*cos(theta/2);
eta = cos((phi+psi)/2)*cos(theta/2);
%%}

%Definition from Wikipedia
%{
eps1 = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
eps2 = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
eps3 = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
eta = cos(phi/2)*cos(phi/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);
%}

%Check that magnitude is 1
quatMag = eps1^2 + eps2^2 + eps3^2 + eta^2;
fprintf('Check Quaternion Magnitude:\n eps1^2 + eps2^2 + eps3^2 + eta^2 =\t %g\n\n',quatMag);

A_quat = zeros(3);

A_quat(1,1) = 1 - 2*(eps2^2 + eps3^2);
A_quat(1,2) = 2*(eps1*eps2 + eps3*eta);
A_quat(1,3) = 2*(eps1*eps3 - eps2*eta);

A_quat(2,1) = 2*(eps2*eps1 - eps3*eta);
A_quat(2,2) = 1 - 2*(eps3^2 + eps1^2);
A_quat(2,3) = 2*(eps2*eps3 + eps1*eta);

A_quat(3,1) = 2*(eps3*eps1 + eps2*eta);
A_quat(3,2) = 2*(eps3*eps2 - eps1*eta);
A_quat(3,3) = 1 - 2*(eps1^2 + eps2^2);


%quiver3(x,y,z,vx,vy,vz,'k');

vx_rot2 = A_quat*vx'; vy_rot2 = A_quat*vy'; vz_rot2 = A_quat*vz';
hold on;
quiver3(x,y,z,vx_rot2',vy_rot2',vz_rot2','b');

axis equal;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
xlabel('x'); ylabel('y'); zlabel('z');
pos = get(gcf, 'Position');
set(gcf,'Position',[600,100,pos(3:4)]);

%}
