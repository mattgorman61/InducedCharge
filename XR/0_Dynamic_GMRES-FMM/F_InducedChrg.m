function [sigma,E_pat,sig_b] = F_InducedChrg(npars,sigma_f_scalar,...
    kappa_air,kappa_p,DeltaArea,NN,x_pat,y_pat,z_pat,NormVec_IF,...
    Ex,Ey,Ez,sig_b)

    lgmres = 1;
    SimpleCount = 0;
    
    NN_tot = NN*npars; % total number of patches
    k_ave = 0.5*(kappa_air+kappa_p);
    k_diff = kappa_air-kappa_p; %kappa_out - kappa_in
    
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
    
    % ES settings: used in FMM
    iprec=2;
    ifcharge=1;
    nsource = NN_tot;
    source = [xc,yc,zc]';
    
    ifdipole=0;
    dipstr = rand(1,nsource);%Not used
    dipvec = rand(3,nsource);%Not used
    ifpot = 0;
    iffld = 1;
    ntarget = 0;
    target = source(:,1:ntarget);%Not used
    target(1,:) = target(1,:) + 10;%Not used
    
    if( ntarget == 0 )
    ifpottarg = 0;
    iffldtarg = 0;
    end
    
    % Iterate to compute induced charge   
    if(lgmres)
    % (1) Compute vector b first
    sig_ds = sig_f.*DeltaArea_Long; % Free charge on each patch: as source in FMM
    charge = sig_ds';
    % call FMM to calculate the field strength by free charge 
    [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
    if( ifpot ), U.pot=U.pot/(4*pi); end
    if( iffld ), U.fld=U.fld/(4*pi); end
    E_f = U.fld';
    % Include external field
    E_f(:,1) = E_f(:,1) + Ex;
    E_f(:,2) = E_f(:,2) + Ey;
    E_f(:,3) = E_f(:,3) + Ez;
    E_f = real(E_f); %keep the real part
    % Projection of E_f on the normal direction
    En_f = E_f(:,1).*NormVec_Long(:,1)+E_f(:,2).*NormVec_Long(:,2)+E_f(:,3).*NormVec_Long(:,3);
    bVec = (1-k_ave)*sig_f-k_diff*En_f;

    % (2)GMRES iteration
    V = zeros(NN_tot); % Store the bases of Krylov subspace
    H = zeros(NN_tot+1,NN_tot); %Hessenburg
    V(:,1) = bVec/norm(bVec);
    
    for i = 1:(NN_tot-1)
        % a...call FMM to compute matrix-vector multiplication
        sig_ds = V(:,i).*DeltaArea_Long; % 
        charge = sig_ds';
        
        [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
        if( ifpot ), U.pot=U.pot/(4*pi); end
        if( iffld ), U.fld=U.fld/(4*pi); end
        
        E_fake = U.fld'; % This is NOT a real field
        E_fake = real(E_fake); %keep the real part
        En_fake = E_fake(:,1).*NormVec_Long(:,1)+E_fake(:,2).*NormVec_Long(:,2)+E_fake(:,3).*NormVec_Long(:,3);
        W=k_ave*V(:,i)+k_diff*En_fake;
        
        % b...Orthogonalization
        for j = 1:i
            H(j,i) = dot(W,V(:,j));
            W = W - H(j,i)*V(:,j);
        end
        % Normalization
        H(j+1,i) = norm(W);
        V(:,i+1) = W/H(j+1,i);
        
        % c...Solve LSQR problem for approaximated solution
        H_sub = H(1:(i+1),1:i);
        eVec = zeros(i+1,1); % eVec = [1,0,...,0]^T with length(eVec) = i+1 in the ith step
        eVec(1) = 1;
        be = norm(bVec)*eVec;
        [y,flag] = lsqr(H_sub,be);
        
        sig_b = V(:,1:i)*y;

        % Check error
        sig_ds = sig_b.*DeltaArea_Long; % 
        charge = sig_ds';
        
        [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
        if( ifpot ), U.pot=U.pot/(4*pi); end
        if( iffld ), U.fld=U.fld/(4*pi); end
    
        E_b = U.fld'; 
        E_b = real(E_b); 
        En_b = E_b(:,1).*NormVec_Long(:,1)+E_b(:,2).*NormVec_Long(:,2)+E_b(:,3).*NormVec_Long(:,3);
        err = bVec-(k_ave*sig_b+k_diff*En_b);
        rel_err = norm(err)/norm(bVec);
%         fprintf('Iteration = %d\n',i);
%         fprintf('Relative error = %f\n',rel_err);
        
        if ((i==(NN_tot-1))&&(rel_err>(1E-4)))
            fprintf('Waring: GMRES does not converge at maximum iteration\n');
            fpringf('Info: i = %d, rel_err = %f\n',i,rel_err);
            pause(10);
        end
        
        if (rel_err<(1E-4))
            fprintf('Iteration = %d\n',i);
            fprintf('Relative error = %f\n',rel_err);
            sigma = sig_b + sig_f;
            E_pat = E_b + E_f;
            break;
        end
        
%         if (max(abs(sig_f))==0) % If neutral
%         sigma = sig_f;
%         E_pat = zeros(NN_tot,3);
%         break;
%         end
    
    end % End of GMRES iteration


    else % if(lgmres)
    % Use simple iteration below
    while(1)
    SimpleCount = SimpleCount + 1;
    sig_b_old = sig_b; % save previous sig_b
    sig = sig_f + sig_b;
    sig_ds = sig.*DeltaArea_Long; % Charge on each patch
    
    % Call FMM to compute E
   
    charge = sig_ds';
    [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
    %[U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
    if( ifpot ), U.pot=U.pot/(4*pi); end
    if( iffld ), U.fld=U.fld/(4*pi); end
    E_str = U.fld';
    % Include external field
    E_str(:,1) = E_str(:,1) + Ex;
    E_str(:,2) = E_str(:,2) + Ey;
    E_str(:,3) = E_str(:,3) + Ez;
    E_str = real(E_str);
%     MaxImag = max(max(imag(E_str)));
%     fprintf('Max imag = %f\n',MaxImag);
    
    % Project E to normal direction
    En = E_str(:,1).*NormVec_Long(:,1)+E_str(:,2).*NormVec_Long(:,2)+E_str(:,3).*NormVec_Long(:,3);
    % Update sig_b
    sig_b = ((1-k_ave)*sig_f - k_diff*En)/k_ave;
    
    % Compute the error
    if (norm(sig_b_old)==0)
        err = 1;
    else
        err = norm(sig_b-sig_b_old)/norm(sig_b_old);
    end
    
    fprintf('Iteration = %d\n',SimpleCount);
    fprintf('Relative error = %f\n',err);
    
    % Check if converge
    if (err < (1E-4))
        sigma = sig_b + sig_f;
    	E_pat = E_str;
        break;
    end
    
%     if (max(sig_f)==0) % If neutral
%         sigma = sig_f;
%         E_pat = zeros(NN_tot,3);
%         break;
%     end
    
    end % while(1) % iteration ends
    
    end % if(lgmres)




