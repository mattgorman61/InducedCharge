function [RMatrx] = F_TransMatrx(QtnEps1,QtnEps2,QtnEps3,QtnEta,npars)
%	Calculate transformation matrix using quaternions
%   Input row vectors: size = i*npars
%   Output Matrix: size = (npars,3,3)
    RMatrx = zeros(npars,3,3);
    
    RMatrx(:,1,1) = 1-2*(QtnEps2.^2 + QtnEps3.^3);
    RMatrx(:,1,2) = 2*(QtnEps1.*QtnEps2 + QtnEps3.*QtnEta);
    RMatrx(:,1,3) = 2*(QtnEps1.*QtnEps3 - QtnEps2.*QtnEta);
    RMatrx(:,2,1) = 2*(QtnEps2.*QtnEps1 - QtnEps3.*QtnEta);
    RMatrx(:,2,2) = 1-2*(QtnEps3.^2 + QtnEps1.^3);
    RMatrx(:,2,3) = 2*(QtnEps2.*QtnEps3 + QtnEps1.*QtnEta);
    RMatrx(:,3,1) = 2*(QtnEps3.*QtnEps1 + QtnEps2.*QtnEta);
    RMatrx(:,3,2) = 2*(QtnEps3.*QtnEps2 - QtnEps1.*QtnEta);
    RMatrx(:,3,3) = 1-2*(QtnEps1.^2 + QtnEps2.^3);
    
end

