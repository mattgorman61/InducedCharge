function [lcollide,Contact1,Contact2,Ematrx_inc] = ...
    F_CollideCheck(i,j,rpar,RMatrx,a,b,c,thickness)
%F_COLLIDECHECK 此处显示有关此函数的摘要
%   此处显示详细说明
    % Translation matrix
    T = zeros(2,4,4);
    T(1,:,:) = [1,0,0,-rpar(i,1); 0,1,0,-rpar(i,2); 0,0,1,-rpar(i,3); 0,0,0,1];
    T(2,:,:) = [1,0,0,-rpar(j,1); 0,1,0,-rpar(j,2); 0,0,1,-rpar(j,3); 0,0,0,1];

    % Rotation matrix
    R = zeros(2,4,4);
    R(1,1:3,1:3) = RMatrx(i,:,:);
    R(1,4,4)=1;
    R(2,1:3,1:3) = RMatrx(j,:,:);
    R(2,4,4)=1;

    % Original ellipsoid matrix
    Eoriginal = diag([1/(a^2),1/(b^2),1/(c^2),-1]); % Original size for contact point
    Eincreased = diag([1/((a+thickness)^2),1/((b+thickness)^2),...
        1/((c+thickness)^2),-1]); % Increased size for collision detection

    % Transform ellipsoid matrix:
    Ematrx_inc = zeros(2,4,4);
    Ematrx_ori = zeros(2,4,4);

    for ii = 1:2
        Tsub = reshape(T(ii,:,:),[4,4]);
        Rsub = reshape(R(ii,:,:),[4,4]);
        Ematrx_inc(ii,:,:)=(Tsub')*(Rsub')...
            *Eincreased*Rsub*Tsub; % Transformed E_increased
        Ematrx_ori(ii,:,:)=(Tsub')*(Rsub')...
            *Eoriginal*Rsub*Tsub;  % Transformed E_original
    end
    
    % Judge if intersect
    Esub1 = reshape(Ematrx_inc(1,:,:),[4,4]);
    Esub2 = reshape(Ematrx_inc(2,:,:),[4,4]);
    Q = inv(Esub1)*Esub2;
    lamda = eig(Q);
    lamda_Img = imag(lamda);
    if(sum(lamda_Img~=0)>0)
        lcollide = 1;
%         fprintf('Collision detected: i=%d, j=%d!\n',i,j);
        
        % Find out collision point here
        [Contact1]=F_CollPoint(i,j,Ematrx_inc);
        [Contact2]=F_CollPoint(j,i,Ematrx_inc);
    else %if(sum(lamda_Img~=0)>0)
        Contact1 = [0,0,0];
        Contact2 = [0,0,0];
        lcollide = 0;
    end %if(sum(lamda_Img~=0)>0)
end

