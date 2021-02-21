function [Fold_par,vpar,vparold,rpar,rot_pf,rotold,thetapar,drotxdtold,...
    drotydtold,drotzdtold] = F_Movement(npars,TransA,dt,mass,...
    Ixpf,Iypf,Izpf,F_par,Fold_par,vpar,vparold,rpar,...
    M_par,rot_pf,rotold,thetapar,drotxdtold,drotydtold,drotzdtold,dampcoeff)
%   Evolve particle movement
%   Translation is evolved in the inertial frame
%   Rotation is evolved in the particle frame
%       

    % (1) Evolve translation
    dvdt = (1.5*F_par-0.5*Fold_par)/mass; % npars*3
    drdt = 1.5*vpar-0.5*vparold;   
    % Update location
    vpar = vpar + dvdt*dt;
    rpar = rpar + drdt*dt;
    % Save info
    Fold_par = F_par;
    vparold = vpar;

    
    % (2) Evolve translation
    % Transform torque
    M_pf = zeros(npars,3);
    for nn = 1:npars
    subA = reshape(TransA(nn,:,:),3,3);
    M_pf = (subA * M_par')';
    end
        
    % Evolve rotation rate
    drotxdt = 1/Ixpf*(rot_pf(:,2).*rot_pf(:,3)*(Iypf-Izpf)+M_pf(:,1));
    drotydt = 1/Iypf*(rot_pf(:,3).*rot_pf(:,1)*(Izpf-Izpf)+M_pf(:,2));
    drotzdt = 1/Izpf*(rot_pf(:,1).*rot_pf(:,2)*(Ixpf-Iypf)+M_pf(:,3));
    dthetadt = 1.5*rot_pf-0.5*rotold;
    
    % Update
    rot_pf(:,1) = rot_pf(:,1) + (1.5*drotxdt-0.5*drotxdtold)*dt;
    rot_pf(:,2) = rot_pf(:,2) + (1.5*drotydt-0.5*drotydtold)*dt;
    rot_pf(:,3) = rot_pf(:,3) + (1.5*drotzdt-0.5*drotzdtold)*dt;
    
    rot_pf(:,1) = (1-dampcoeff)*rot_pf(:,1);
    rot_pf(:,2) = (1-dampcoeff)*rot_pf(:,2);
    rot_pf(:,3) = (1-dampcoeff)*rot_pf(:,3);
    
    thetapar = thetapar + dthetadt*dt;
    
    % Save old info
    drotxdtold = drotxdt;
    drotydtold = drotydt;
    drotzdtold = drotzdt;
    rotold = rot_pf;
        
end

