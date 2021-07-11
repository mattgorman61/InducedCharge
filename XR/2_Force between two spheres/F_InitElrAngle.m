function [Elr_phi,Elr_psi,Elr_theta] = F_InitElrAngle(thetax,thetay,thetaz)
% Convert initial orientation angle to Euler angle
%   Follow the definition in DEM
    Elr_phi = 0.5*thetaz;
    Elr_psi = 0.5*thetay;
    Elr_theta = 0.5*thetax;
end

