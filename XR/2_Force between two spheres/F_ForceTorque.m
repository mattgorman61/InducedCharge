function [F_par,M_par] = F_ForceTorque(sigma,E_pat,DeltaArea,NN,npars,...
    x_pat,y_pat,z_pat,xpar,ypar,zpar)
%F_FORCETORQUE 此处显示有关此函数的摘要
%   此处显示详细说明
    % Definition
    F_pat = zeros(NN,npars,3); % (patch id,particle id,dimension)
    M_pat = zeros(NN,npars,3); % (patch id,particle id,dimension)
    F_par = zeros(npars,3);
    M_par = zeros(npars,3);
    
    % Relative location of each patch on each particle
    xx_pat = reshape(x_pat,[NN,npars]);
    yy_pat = reshape(y_pat,[NN,npars]);
    zz_pat = reshape(z_pat,[NN,npars]);
    
    xx_pat_rel = zeros(NN,npars);
    yy_pat_rel = zeros(NN,npars);
    zz_pat_rel = zeros(NN,npars);
    
    for ii = 1:npars
    xx_pat_rel(:,ii) = xx_pat(:,ii)-xpar(ii);
    yy_pat_rel(:,ii) = yy_pat(:,ii)-ypar(ii);
    zz_pat_rel(:,ii) = zz_pat(:,ii)-zpar(ii);
    end
    
    x_rel_pat = reshape(xx_pat_rel,[],1);
    y_rel_pat = reshape(yy_pat_rel,[],1);
    z_rel_pat = reshape(zz_pat_rel,[],1);
    
    DeltaArea_Long = zeros(NN*npars,1);
    for i = 1:npars
        DeltaArea_Long(((i-1)*NN+1):(i*NN)) = DeltaArea;
    end
    
    % Force on each patch
    FFx_pat = E_pat(:,1).*sigma.*DeltaArea_Long;
    FFy_pat = E_pat(:,2).*sigma.*DeltaArea_Long;
    FFz_pat = E_pat(:,3).*sigma.*DeltaArea_Long;
    
    % Torque on each patch
    r_rel_pat = [x_rel_pat,y_rel_pat,z_rel_pat];
    FF_pat = [FFx_pat,FFy_pat,FFz_pat];
    MM_pat = cross(r_rel_pat,FF_pat,2);
    
    % Reshape force on each patch
    F_pat(:,:,1) = reshape(FFx_pat,[NN,npars]);
    F_pat(:,:,2) = reshape(FFy_pat,[NN,npars]);
    F_pat(:,:,3) = reshape(FFz_pat,[NN,npars]);
    
    % Reshape torque on each patch
    M_pat(:,:,1) = reshape(MM_pat(:,1),[NN,npars]);
    M_pat(:,:,2) = reshape(MM_pat(:,2),[NN,npars]);
    M_pat(:,:,3) = reshape(MM_pat(:,3),[NN,npars]);
       
    % Force/Torque on each particle
    for ii = 1:npars
    F_par(ii,1) = sum(F_pat(:,ii,1));
    F_par(ii,2) = sum(F_pat(:,ii,2));
    F_par(ii,3) = sum(F_pat(:,ii,3));
    
    M_par(ii,1) = sum(M_pat(:,ii,1));
    M_par(ii,2) = sum(M_pat(:,ii,2));
    M_par(ii,3) = sum(M_pat(:,ii,3));
    end
    
end

