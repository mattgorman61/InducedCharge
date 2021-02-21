function [ContactPoint] = F_CollPoint(mm,nn,Ematrx)
% Find out the contact point in case of collisions
%  mm = id of the targer particle (Matrix: q_ij)
%  nn = id of the other particle (Matrix: p_ij)
    
% Characteristic matrix of ellipsoid
q = reshape(Ematrx(mm,:,:),[4,4]); % target particle
p = reshape(Ematrx(nn,:,:),[4,4]); % the other one

% q and p are not strictly symetric because of numerical errors
% So I pick out the symmetric part
% The symmetric part is close to its original form
q = 0.5*(q+q');
p = 0.5*(p+p');

% Calculate the 6-order polynomial of tau
a1 = p(2,2)*p(3,3)-p(2,3)^2;
a2 = q(2,2)*p(3,3)+p(2,2)*q(3,3)-2*p(2,3)*q(2,3);
a3 = q(2,2)*q(3,3)-q(2,3)^2;
a4 = p(1,2)*p(3,3)-p(1,3)*p(2,3);
a5 = q(1,2)*p(3,3)+p(1,2)*q(3,3)-p(1,3)*q(2,3)-q(1,3)*p(2,3);
a6 = q(1,2)*q(3,3)-q(1,3)*q(2,3);
a7 = p(1,2)*p(2,3)-p(1,3)*p(2,2);
a8 = q(1,2)*p(2,3)+p(1,2)*q(2,3)-p(1,3)*q(2,2)-q(1,3)*p(2,2);
a9 = q(1,2)*q(2,3)-q(1,3)*q(2,2);

b1 = q(1,1)*a1+p(1,1)*a2;
b2 = p(1,1)*a3+q(1,1)*a2;
b3 = q(1,1)*a3;
b4 = q(1,2)*a4+p(1,2)*a5;
b5 = p(1,2)*a6+q(1,2)*a5;
b6 = q(1,2)*a6;
b7 = q(1,3)*a7+p(1,3)*a8;
b8 = q(1,3)*a8+p(1,3)*a9;
b9 = q(1,3)*a9;

c1 = p(1,1)*a1-p(1,2)*a4+p(1,3)*a7;
c2 = b1-b4+b7;
c3 = b2-b5+b8;
c4 = b3-b6+b9;

d1 = 2*c1*c2;
d2 = 2*c1*c3+c2^2;
d3 = 2*(c1*c4+c2*c3);
d4 = 2*c2*c4+c3^2;
d5 = 2*c3*c4;

aa1 = p(2,2)*p(3,3)-p(2,3)^2;
aa2 = p(2,2)*q(3,3)+q(2,2)*p(3,3)-2*p(2,3)*q(2,3);
aa3 = q(2,2)*q(3,3)-q(2,3)^2;
aa4 = p(1,3)*p(2,3)-p(1,2)*p(3,3);
aa5 = p(2,3)*q(1,3)+p(1,3)*q(2,3)-p(3,3)*q(1,2)-p(1,2)*q(3,3);
aa6 = q(1,3)*q(2,3)-q(1,2)*q(3,3);
aa7 = p(1,2)*p(2,3)-p(2,2)*p(1,3);
aa8 = p(1,2)*q(2,3)+p(2,3)*q(1,2)-p(1,3)*q(2,2)-p(2,2)*q(1,3);
aa9 = q(1,2)*q(2,3)-q(2,2)*q(1,3);
aa10 = p(1,1)*p(3,3)-p(1,3)^2;
aa11 = p(3,3)*q(1,1)+p(1,1)*q(3,3)-2*p(1,3)*q(1,3);
aa12 = q(1,1)*q(3,3)-q(1,3)^2;
aa13 = p(1,2)*p(1,3)-p(1,1)*p(2,3);
aa14 = p(1,2)*q(1,3)+p(1,3)*q(1,2)-p(1,1)*q(2,3)-p(2,3)*q(1,1);
aa15 = q(1,2)*q(1,3)-q(1,1)*q(2,3);
aa16 = p(1,1)*p(2,2)-p(1,2)^2;
aa17 = p(1,1)*q(2,2)+p(2,2)*q(1,1)-2*p(1,2)*q(1,2);
aa18 = q(1,1)*q(2,2)-q(1,2)^2;

bb1 = p(1,4)*aa1+p(2,4)*aa4+p(3,4)*aa7;
bb2 = p(1,4)*aa2+p(2,4)*aa5+p(3,4)*aa8 + ...
      q(1,4)*aa1+q(2,4)*aa4+q(3,4)*aa7;
bb3 = p(1,4)*aa3+p(2,4)*aa6+p(3,4)*aa9 + ...
      q(1,4)*aa2+q(2,4)*aa5+q(3,4)*aa8;
bb4 = q(1,4)*aa3+q(2,4)*aa6+q(3,4)*aa9;
bb5 = p(1,4)*aa4+p(2,4)*aa10+p(3,4)*aa13;
bb6 = p(1,4)*aa5+p(2,4)*aa11+p(3,4)*aa14 + ...
      q(1,4)*aa4+q(2,4)*aa10+q(3,4)*aa13;
bb7 = p(1,4)*aa6+p(2,4)*aa12+p(3,4)*aa15 + ...
      q(1,4)*aa5+q(2,4)*aa11+q(3,4)*aa14;
bb8 = q(1,4)*aa6+q(2,4)*aa12+q(3,4)*aa15;
bb9 = p(1,4)*aa7+p(2,4)*aa13+p(3,4)*aa16;
bb10 = p(1,4)*aa8+p(2,4)*aa14+p(3,4)*aa17 + ...
       q(1,4)*aa7+q(2,4)*aa13+q(3,4)*aa16;
bb11 = p(1,4)*aa9+p(2,4)*aa15+p(3,4)*aa18 + ...
       q(1,4)*aa8+q(2,4)*aa14+q(3,4)*aa17;
bb12 = q(1,4)*aa9+q(2,4)*aa15+q(3,4)*aa18;

cc1 = q(1,4)*bb1+q(2,4)*bb5+q(3,4)*bb9;
cc2 = q(1,4)*bb2+q(2,4)*bb6+q(3,4)*bb10;
cc3 = q(1,4)*bb3+q(2,4)*bb7+q(3,4)*bb11;
cc4 = q(1,4)*bb4+q(2,4)*bb8+q(3,4)*bb12;

dd1 = c1*cc1;
dd2 = c1*cc2+c2*cc1;
dd3 = c1*cc3+c2*cc2+c3*cc1;
dd4 = c1*cc4+c2*cc3+c3*cc2+c4*cc1;
dd5 = c2*cc4+c3*cc3+c4*cc2;
dd6 = c3*cc4+c4*cc3;
dd7 = c4*cc4;

ee1 = 2*bb2*bb1;
ee2 = 2*bb1*bb3+bb2^2;
ee3 = 2*bb1*bb4+2*bb2*bb3;
ee4 = 2*bb2*b4+bb3^2;
ee5 = 2*bb3*bb4;
ee6 = 2*bb5*bb6;
ee7 = 2*bb5*bb7+bb6^2;
ee8 = 2*bb5*bb8+2*bb6*bb7;
ee9 = 2*bb6*bb8+bb7^2;
ee10 = 2*bb7*bb8;
ee11 = 2*bb9*bb10;
ee12 = 2*bb9*bb11+bb10^2;
ee13 = 2*bb9*bb12+2*bb10*bb11;
ee14 = 2*bb10*bb12+bb11^2;
ee15 = 2*bb11*bb12;
ee16 = bb1*bb6+bb2*bb5;
ee17 = bb1*bb7+bb5*bb3+bb2*bb6;
ee18 = bb1*bb8+bb4*bb5+bb2*bb7+bb3*bb6;
ee19 = bb2*bb8+bb4*bb6+bb3*bb7;
ee20 = bb3*bb8+bb4*bb7;
ee21 = bb1*bb10+bb2*bb9;
ee22 = bb1*bb11+bb3*bb9+bb2*bb10;
ee23 = bb1*bb12+bb9*bb4+bb2*bb11+bb3*bb10;
ee24 = bb2*bb12+bb4*bb10+bb3*bb11;
ee25 = bb3*bb12+bb4*bb11;
ee26 = bb5*bb10+bb6*bb9;
ee27 = bb5*bb11+bb9*bb7+bb6*bb10;
ee28 = bb5*bb12+bb6*bb11+bb7*bb10+bb8*bb9;
ee29 = bb6*bb12+bb8*bb10+bb7*bb11;
ee30 = bb7*bb12+bb8*bb11;

ff1 = q(1,1)*bb1^2+q(2,2)*bb5^2+q(3,3)*bb9^2+2*q(1,2)*bb1*bb5+...
      2*q(1,3)*bb1*bb9+2*q(2,3)*bb5*bb9;
ff2 = q(1,1)*ee1+q(2,2)*ee6+q(3,3)*ee11+2*q(1,2)*ee16+2*q(1,3)*ee21+2*q(2,3)*ee26;
ff3 = q(1,1)*ee2+q(2,2)*ee7+q(3,3)*ee12+2*q(1,2)*ee17+2*q(1,3)*ee22+2*q(2,3)*ee27;
ff4 = q(1,1)*ee3+q(2,2)*ee8+q(3,3)*ee13+2*q(1,2)*ee18+2*q(1,3)*ee23+2*q(2,3)*ee28;
ff5 = q(1,1)*ee4+q(2,2)*ee9+q(3,3)*ee14+2*q(1,2)*ee19+2*q(1,3)*ee24+2*q(2,3)*ee29;
ff6 = q(1,1)*ee5+q(2,2)*ee10+q(3,3)*ee15+2*q(1,2)*ee20+2*q(1,3)*ee25+2*q(2,3)*ee30;
ff7 = q(1,1)*bb4^2+q(2,2)*bb8^2+q(3,3)*bb12^2+2*q(1,2)*bb4*bb8+2*q(1,3)*bb4*bb12...
    +2*q(2,3)*bb8*bb12;

gg7 = ff7-2*dd7+q(4,4)*c4^2;
gg6 = ff6-2*dd6+q(4,4)*d5;
gg5 = ff5-2*dd5+q(4,4)*d4;
gg4 = ff4-2*dd4+q(4,4)*d3;
gg3 = ff3-2*dd3+q(4,4)*d2;
gg2 = ff2-2*dd2+q(4,4)*d1;
gg1 = ff1-2*dd1+q(4,4)*c1^2;

% Find the roots of the following polynomial:
% gg7*tau^6 + gg6*tau^5 + gg5*tau^4 + gg4*tau^3 + gg3*tau^2 + gg2*tau^1 ...
% + gg1 = 0

TauPolyNom = [gg7,gg6,gg5,gg4,gg3,gg2,gg1];

RootTau = roots(TauPolyNom);

A1 = Ematrx(mm,1:3,1:3);
A2 = Ematrx(nn,1:3,1:3);
bVec1 = 2*Ematrx(mm,1:3,4);
bVec2 = 2*Ematrx(nn,1:3,4);

% Calculate and check the contact point
alpha_min = 999; % Large enough

for ii = 1:length(RootTau)
    if (imag(RootTau(ii))==0)
        AMatrx = A2+RootTau(ii)*A1;
        AMatrx = reshape(AMatrx,[3,3]);
        bVec = (bVec2+RootTau(ii)*bVec1)';
        Point = -0.5*inv(AMatrx)*bVec; % Possible contact point
        % Find smallest alpha and corresponding contact point
        alpha = Point'*reshape(A2,[3,3])*Point+bVec2*Point+q(4,4);
        if (alpha<alpha_min)
            alpha_min = alpha;
            ContactPoint = Point;
        end
    end %if (imag(RootTau(ii))==0)
end % for ii = 1:length(RootTau)

fprintf('Contact point on id=%d is: x=%f,y=%f,z=%f\n',mm,...
    ContactPoint(1),ContactPoint(2),ContactPoint(3));
    
    ContactPoint = ContactPoint';
    
end

