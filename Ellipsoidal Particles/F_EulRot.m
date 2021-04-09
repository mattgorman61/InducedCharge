function [x,y,z,A_Eul] = F_EulRot(x0,y0,z0,psi,phi,theta)
% Rotate Ellipsoid using Euler Angles
%{
Given:
    
    x0.................... original vector of x locations
    y0.................... original vector of y locations
    z0.................... original vector of z locations
    psi................... Euler rotation angle
    phi................... Euler rotation angle
    theta................. Euler rotation angle
    
   
Returns:
    x..................... rotated vector of x-locations
    y..................... rotated vector of y-locations
    z..................... rotated vector of z-locations
    dA.................... vector of Area elements for each patch
    A_Eul................. rotation Matrix



%}

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


posMat = [x0,y0,z0]';
rot_posMat = A_Eul*posMat;
rot_posMat = rot_posMat';

x = rot_posMat(:,1); y = rot_posMat(:,2); z = rot_posMat(:,3);




end

