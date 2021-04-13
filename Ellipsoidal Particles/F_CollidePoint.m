function [x_coll,y_coll,z_coll] = F_CollidePoint(a0,b0,c0,x0,y0,z0,AEul_0, a1,b1,c1,x1,y1,z1,AEul_1)
% CHECKS WHETHER TWO ELLIPSES HAVE COLLIDED
%   Alfano and Greer 2003 Journ. of Guidance, Control, and Dynamics


lCollideFlag = -1;

% First Ellipse
A0 = 1/a0^2; B0 = 1/b0^2; C0 = 1/c0^2; 
D0 = 0; E0 = 0; F0 = 0;
G0 = -2*x0/a0^2; H0 = -2*y0/b0^2; J0 = -2*z0/c0^2;
K0 = (x0^2)/a0^2 + (y0^2)/b0^2 + (z0^2)/c0^2 - 1;

S_0 = [2*A0, D0, F0, G0; D0, 2*B0, E0, H0; F0, E0, 2*C0, J0; G0,H0,J0,2*K0];
S_0 = AEul_0*S_0*AEul_0'; % Apply Rotation

% Second Ellipse
A1 = 1/a1^2; B1 = 1/b1^2; C1 = 1/c1^2; 
D1 = 0; E1 = 0; F1 = 0;
G1 = -2*x1/a1^2; H1 = -2*y1/b1^2; J1 = -2*z1/c1^2;
K1 = (x1^2)/a1^2 + (y1^2)/b1^2 + (z1^2)/c1^2 - 1;

S_1 = [2*A1, D1, F1, G1; D1, 2*B1, E1, H1; F1, E1, 2*C1, J1; G1,H1,J1,2*K1];
S_1 = AEul_1*S_1*AEul_1'; % Apply Rotation

EigMat = (S_0)\S_1;
[EigVect,EigValues] = eig(EigMat);
[~,cnum] = size(EigVect);

for i=1:cnum
    evecti = EigVect(:,i);
    evali = EigValues(:,i);
    if(isequal(evecti,[x0;y0;z0;1]))
        lCollideFlag = 1;
    end
    if(~isequal(imag(evecti),[0;0;0;0]))
        lCollideFlag = 1;
    end
    if(~isequal(imag(evali),[0;0;0;0]))
        lCollideFlag = 1;
    end
end

if(lCollideFlag > 0)
    lCollide = true;
    
else
    lCollide = false;
end

end

