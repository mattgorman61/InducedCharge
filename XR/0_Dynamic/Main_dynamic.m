clear;
clc;

% =========================================================================
% ========== Definition ==========
lCheckGeometry = 0;
lCheckPotential = 0;
dampcoeff = 0.0;

filename = 'a2b1c1.mat'; % Geometry: Axisymmetric ellipsoid

kappa_air = 1; % Dielectric constance of air
kappa_p = 2.5; % Dielectric constance of particle
epsilon = 1; %Unit: F/m
rhop = 1; % particle density

dt = 0.002;
itmax = 1000;
thickness = 0.0; % Preset gap in collision detection

MyThetaCheck = zeros(itmax,3);

%==========================================================================
% Initialization
% Vectors are used to store info of all particles
% Note: vector length should equal npars
npars = 2;

sigma_f_scalar = [1,-1];%linspace(0,0,npars); % Surface free charge density
Ex = 0; % External field
Ey = 0;
Ez = 0;

% Location
xpar = [0,4];%linspace(0,0,npars); % size=1*npars
ypar = [0,2];
zpar = [0,2];%linspace(0,0,npars);
rpar = [xpar;ypar;zpar]'; % size=npars*3

% Velocity
vpar = zeros(npars,3); % size=npars*3
vparold = zeros(npars,3); % size=npars*3

% Orientation angle
thetapar = zeros(npars,3);

% Rotation rate
rot = zeros(npars,3);
rotold = zeros(npars,3);
drotxdtold = zeros(npars,1);
drotydtold = zeros(npars,1);
drotzdtold = zeros(npars,1);

% Force/Torque
F_par = zeros(npars,3);
Fold_par = zeros(npars,3);
M_par = zeros(npars,3);
Mold_par = zeros(npars,3);
M_pf = zeros(npars,3); % particle frame
drotxdt = zeros(npars,3); % particle frame
drotydt = zeros(npars,3); % particle frame
drotzdt = zeros(npars,3); % particle frame

% Quaternion rate
Eps1rateold = zeros(npars,1);
Eps2rateold = zeros(npars,1);
Eps3rateold = zeros(npars,1);
Etarateold = zeros(npars,1);

% Read patch info
[x_rel,y_rel,z_rel,DeltaArea,NN,NormVec,a,b,c]=F_GeometryEllip(filename);

mass = 4/3*pi*a*b*c; % particle mass
Ixpf = mass*(b^2+c^2)/5; % moment of inertia in particle frame
Iypf = mass*(a^2+c^2)/5; % moment of inertia in particle frame
Izpf = mass*(a^2+b^2)/5; % moment of inertia in particle frame

% Initial Euler angle
[Elr_phi,Elr_psi,Elr_theta] = F_InitElrAngle(thetapar(:,1),thetapar(:,2),...
    thetapar(:,3));

% Euler angles to quaternions
[QtnEps1,QtnEps2,QtnEps3,QtnEta] = F_Elr2Qtn(Elr_phi,Elr_psi,Elr_theta);

% Initial transformation matrix from quaternions
[RMatrx] = F_TransMatrx(QtnEps1,QtnEps2,QtnEps3,QtnEta,npars);

% Transform rotation rate
rot_pf = zeros(npars,3);
for nn = 1:npars
    subA = reshape(RMatrx(nn,:,:),3,3);
    rot_pf = (subA * rot')';
end

% Temporal evolution starts
for tt = 1:itmax
    
    fprintf('Time step tt = %d\n',tt);
    
    % (1) Update patch location in the inertial frame
    
    [x_pat,y_pat,z_pat,NormVec_IF] = F_AbsPatch(xpar,ypar,zpar,x_rel,y_rel,...
    z_rel,RMatrx,NormVec,npars,NN);
    
    % (2) Calculate induced surface charge
    [sigma,E_pat]=F_InducedChrg(npars,sigma_f_scalar,...
        kappa_air,kappa_p,DeltaArea,NN,x_pat,y_pat,z_pat,NormVec_IF,...
        Ex,Ey,Ez);
    
    % (3) Calculate force and torque
    % Note: F_pat(patch id,particle id,dimension)
    [F_par,M_par] = F_ForceTorque(sigma,E_pat,DeltaArea,NN,npars,...
    x_pat,y_pat,z_pat,xpar,ypar,zpar);
% 
%     fprintf('Fx_par = %f, Mx_par = %f\n',F_par(1,1),M_par(1,1));
%     fprintf('Fy_par = %f, My_par = %f\n',F_par(1,2),M_par(1,2));
%     fprintf('Fz_par = %f, Mz_par = %f\n',F_par(1,3),M_par(1,3));

    % (4) Evolve ellipsoid movements
    
    % Update particle movements
    [Fold_par,vpar,vparold,rpar,rot_pf,rotold,thetapar,drotxdtold,...
    drotydtold,drotzdtold] = F_Movement(npars,RMatrx,dt,mass,...
    Ixpf,Iypf,Izpf,F_par,Fold_par,vpar,vparold,rpar,...
    M_par,rot_pf,rotold,thetapar,drotxdtold,drotydtold,drotzdtold,dampcoeff);
    
    % Update quaternion
    [QtnEps1,QtnEps2,QtnEps3,QtnEta, Eps1rateold,Eps2rateold,...
    Eps3rateold,Etarateold] = F_EvlvQtn(dt,QtnEps1,QtnEps2,QtnEps3,QtnEta,...
    Eps1rateold,Eps2rateold,Eps3rateold,Etarateold,rot_pf);
    
    % Update transformation matrix
    [RMatrx] = F_TransMatrx(QtnEps1,QtnEps2,QtnEps3,QtnEta,npars);
    
    % Update location
    xpar = rpar(:,1)'; % size=1*npars
    ypar = rpar(:,2)';
    zpar = rpar(:,3)';
            
    % Plot charge distribution of present time step
    figure(99);
    
    xplot = reshape(x_pat,[npars*NN,1]);
    yplot = reshape(y_pat,[npars*NN,1]);
    zplot = reshape(z_pat,[npars*NN,1]);
    scatter3(xplot,yplot,zplot,10,sigma,'filled');
        
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(gca,'LineWidth',1.5);
    axis equal;
    
    % Collision detection in multiple-particle case
    if (npars >1)
        % Construct collision detection list
        [DetectList] = F_DetectList(npars,rpar,a,b,c);
        
        % Check possible collision pairs & find out collision point
        [CollideList,ContactPoints] = F_Collision(npars,DetectList,a,b,c,...
    rpar,RMatrx,thickness);
        
        if (sum(CollideList)>0)
            fprintf('Collision detected!\n');
            % Plot snapshot
            [x_pat,y_pat,z_pat,NormVec_IF] = F_AbsPatch(xpar,ypar,zpar,x_rel,y_rel,...
    z_rel,RMatrx,NormVec,npars,NN);
            figure(98);
            xplot = reshape(x_pat,[npars*NN,1]);
            yplot = reshape(y_pat,[npars*NN,1]);
            zplot = reshape(z_pat,[npars*NN,1]);
            scatter3(xplot,yplot,zplot,10,sigma,'filled');
            hold on;
            
            % Plot the contact point
            scatter3(ContactPoints(1,2,1),ContactPoints(1,2,2),...
            ContactPoints(1,2,3),20,'r','filled');
            scatter3(ContactPoints(2,1,1),ContactPoints(2,1,2),...
            ContactPoints(2,1,3),20,'b','filled');
            
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            set(gca,'LineWidth',1.5);
            axis equal;
            
            break;
        end % if (sum(CollideList)>0)

        % Update post-collision velocity/rotation rate
        % Based on hard sphere model
        
    end % if (npars >1)

end % for tt = 1:itmax 
