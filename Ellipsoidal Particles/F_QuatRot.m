function [x,y,z] = F_EulRot(x0,y0,z0,psi,phi,theta)
% Rotate Ellipsoid using Quaternions

eps1 = cos((phi-psi)/2)*sin(theta/2); 
eps2 = sin((phi-psi)/2)*sin(theta/2);
eps3 = sin((phi+psi)/2)*cos(theta/2);
eta = cos((phi+psi)/2)*cos(theta/2);

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


posMat = [x0,y0,z0]';
rot_posMat = A_quat*posMat;
rot_posMat = rot_posMat';

x = rot_posMat(:,1); y = rot_posMat(:,2); z = rot_posMat(:,3);




end

