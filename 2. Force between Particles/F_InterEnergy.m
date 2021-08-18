function [E_self] = F_InterEnergy(npars,NN,sigma_f_scalar,...
    DeltaArea,sigma,x_pat,y_pat,z_pat,Correct)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

    sig_f = Correct*sigma_f_scalar .* linspace(1,1,NN)';
    
    % Compute self energy of each particle
    % ES settings: used in FMM
    iprec=2;
    ifcharge=1;
    nsource = NN;
        
    ifdipole=0;
    dipstr = rand(1,nsource);%Not used
    dipvec = rand(3,nsource);%Not used
    ifpot = 1;
    iffld = 0;
    ntarget = 0;
        
    if( ntarget == 0 )
    ifpottarg = 0;
    iffldtarg = 0;
    end
    
    sigma_par = reshape(sigma,[NN,npars]);
    E_self = zeros(npars,1)';
    
    for ii = 1:npars
        % Source info
        source = [x_pat(:,ii),y_pat(:,ii),z_pat(:,ii)]';
        target = source(:,1:ntarget);%Not used
        target(1,:) = target(1,:) + 10;%Not used
        sig_ds = sigma_par(:,ii).*DeltaArea;
        charge = sig_ds';
        % Call FMM
        [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
        if( ifpot ), U.pot=U.pot/(4*pi); end
        if( iffld ), U.fld=U.fld/(4*pi); end
        
        Psi = real(U.pot)';
        E_self(ii) = 0.5*sum(Psi.*sig_f(:,ii).*DeltaArea);
    end
end

