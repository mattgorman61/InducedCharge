function [sigma,E_pat] = F_InducedChrg(npars,sigma_f_scalar,...
    kappa_air,kappa_p,DeltaArea,NN,x_pat,y_pat,z_pat,NormVec_IF,...
    Ex,Ey,Ez)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    
    NN_tot = NN*npars; % total number of patches
    
    xc = reshape(x_pat,[],1);
    yc = reshape(y_pat,[],1);
    zc = reshape(z_pat,[],1);
    
    sig_f = reshape(sigma_f_scalar .* linspace(1,1,NN)',[],1);
    
    NormVec_Long = zeros(NN_tot,3);
    DeltaArea_Long = zeros(NN_tot,1);
    for i = 1:npars
        NormVec_Long(((i-1)*NN+1):(i*NN),:) = NormVec_IF(i,:,:);
        DeltaArea_Long(((i-1)*NN+1):(i*NN)) = DeltaArea;
    end
    
%     rdiff = zeros(NN,NN,3);
    rdiff(:,:,1) = xc - xc';
    rdiff(:,:,2) = yc - yc';
    rdiff(:,:,3) = zc - zc';

    rdist1 = pdist([xc,yc,zc]);
    rdist = squareform(rdist1);

    I = rdiff ./(rdist.^3)/4/pi; % Diagonal components become NaN
    for i = 1:NN_tot
        I(i,i,:) = [0,0,0];
    end

    L = zeros(NN_tot,NN_tot);
    for i = 1:NN_tot
        L(i,:) = I(i,:,1)*NormVec_Long(i,1)+I(i,:,2)*NormVec_Long(i,2)+...
            I(i,:,3)*NormVec_Long(i,3);
    end

    k_ave = 0.5*(kappa_air+kappa_p);
    k_diff = kappa_air-kappa_p; %kappa_out - kappa_in

    dSj = DeltaArea_Long.*linspace(1,1,NN_tot)';% Assume patch area is the same,should be updated if necessary
    A = k_ave*eye(NN_tot)+k_diff*L*dSj(1); % Left matrix
    B = (1-k_ave)*eye(NN_tot)-k_diff*L*dSj(1); % Geometry matrix for free charge

    % Compute contribution by external field
    E_proj=linspace(0,0,NN_tot)'; % No external field
    for i = 1:NN_tot
    E_proj(i)=Ex*NormVec_Long(i,1)+Ey*NormVec_Long(i,2)+Ez*NormVec_Long(i,3);
    end

    b = B*sig_f - k_diff*E_proj;
    tolerance = norm(b,2)/10000;
    % use GMRES to iterate
    sig_b = gmres(A,b,10,tolerance);
    % total charge density
    sigma = sig_f + sig_b;
    
    % Calculate field strength on each patch
    E_pat = zeros(NN_tot,3);
    
    E_pat(:,1) = I(:,:,1)*sigma + Ex;
    E_pat(:,2) = I(:,:,2)*sigma + Ey;
    E_pat(:,3) = I(:,:,3)*sigma + Ez;
    
    % Display the surface charge distribution
    
%     figure(3);
%     size = 10;
%     scatter3(xc(:),yc(:),zc(:),size,sigma(:),'filled');
%     hold on;
%     axis equal;
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     set(gca,'LineWidth',1.5);

end

