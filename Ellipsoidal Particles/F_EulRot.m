function [x,y,z] = F_EulRot(x0,y0,z0,psi,phi,theta)
% Rotate Ellipsoid using Euler Angles

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


posMat = [x0,y0,z0]';
rot_posMat = A_Eul*posMat;
rot_posMat = rot_posMat';

x = rot_posMat(:,1); y = rot_posMat(:,2); z = rot_posMat(:,3);




end

