clc; clear all; close all;


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

%%{
psi = pi/4;
phi = pi/4;
theta = pi/4;

A_Eul = zeros(3);

A_Eul(1,1) = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi);
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

%A_Eul*A_Eul'

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

% A_quat*A_quat'
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
