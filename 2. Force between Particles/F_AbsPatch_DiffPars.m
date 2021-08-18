function [x_pat,y_pat,z_pat,NormVec_IF] = F_AbsPatch_DiffPars(xpar,ypar,zpar,x_rel,y_rel,...
    z_rel,TransA,NormVec,npars,NN)
%   Calculate patch location in the inertial frame
%   Input: xpar,ypar,zpar: 1*npars
%   Input: x_rel,y_rel,z_rel: npatch*1
%   Input: TransA: npars*3*3
%   Input: NormVec: npatch*3
%   Output: x_pat,y_pat,z_pat: npatch*npars
%   Output: NormVec_IF: TBD
    
    x_pat = zeros(NN,npars);
    y_pat = zeros(NN,npars);
    z_pat = zeros(NN,npars);
    NormVec_IF = zeros(npars,NN,3);
    
    % Note: from particle frame to inertial frame, the trans matrix is the
    % transpose of A (i.e. A')
    for i = 1:npars
            subA = reshape(TransA(i,:,:),3,3);
            r_rel=(subA'*[x_rel,y_rel,z_rel]')';
            x_pat(:,i)=xpar(i)*linspace(1,1,NN)'+r_rel(:,1);
            y_pat(:,i)=ypar(i)*linspace(1,1,NN)'+r_rel(:,2);
            z_pat(:,i)=zpar(i)*linspace(1,1,NN)'+r_rel(:,3);
            NormVec_IF(i,:,:)=(subA'*(NormVec'))';
    end    
end

