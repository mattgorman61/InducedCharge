function [QtnEps1,QtnEps2,QtnEps3,QtnEta] = F_Elr2Qtn(Elr_Phi,Elr_Psi,Elr_Theta)
%   Convert Euler angles to quaternions
%   Input row vectors: size = 1*npars
%   Output row vectors: size = 1*npars
    
%     If follow definition in paper, Eqs. should be as follows
%     QtnEps1 = cos(0.5*(Elr_Phi-Elr_Psi)).*sin(0.5*Elr_Theta);
%     QtnEps2 = sin(0.5*(Elr_Phi-Elr_Psi)).*sin(0.5*Elr_Theta);
%     QtnEps3 = sin(0.5*(Elr_Phi+Elr_Psi)).*cos(0.5*Elr_Theta);
%     QtnEta = cos(0.5*(Elr_Phi+Elr_Psi)).*cos(0.5*Elr_Theta);

%     Follow the definition in DEM, Eqs. become as follows
    QtnEps1 = cos(Elr_Phi-Elr_Psi).*sin(Elr_Theta);
    QtnEps2 = sin(Elr_Phi-Elr_Psi).*sin(Elr_Theta);
    QtnEps3 = sin(Elr_Phi+Elr_Psi).*cos(Elr_Theta);
    QtnEta = cos(Elr_Phi+Elr_Psi).*cos(Elr_Theta);

end

