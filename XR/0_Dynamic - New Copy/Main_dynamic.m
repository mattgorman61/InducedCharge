clear;
clc;

% =========================================================================
% ========== Definition ==========
lCheckGeometry = 0;
lCheckPotential = 0;
erest = 0.9;
frict = 0.3;
Elastic = 1E5;
Poisson = 0.33;

filename = 'a1b1c1.mat'; % Geometry: Axisymmetric ellipsoid

kappa_air = 1; % Dielectric constance of air
kappa_p = 2.5; % Dielectric constance of particle
epsilon = 1; %Unit: F/m
rhop = 1; % particle density

dt = 1E-3; % 1E-3
itmax = 3000; % 10000
thickness = 0.2; % Preset gap in collision detection

figure_id = 1;

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
xpar = [-1,2];%linspace(0,0,npars); % size=1*npars
ypar = [0,0];
zpar = [0,0];%linspace(0,0,npars);
rpar = [xpar;ypar;zpar]'; % size=npars*3

% Velocity
vpar = [0,0,0;-2,2,0]; % size=npars*3
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
I_pf = diag([Ixpf,Iypf,Izpf]); % inertia tensor in particle frame

% Initial Euler angle
[Elr_phi,Elr_psi,Elr_theta] = F_InitElrAngle(thetapar(:,1),thetapar(:,2),...
    thetapar(:,3));

% Euler angles to quaternions
[QtnEps1,QtnEps2,QtnEps3,QtnEta] = F_Elr2Qtn(Elr_phi,Elr_psi,Elr_theta);

% Initial transformation matrix from quaternions
[RMatrx] = F_TransMatrx(QtnEps1,QtnEps2,QtnEps3,QtnEta,npars);

% Transform rotation rate
rot_pf = zeros(npars,3);
rotold_pf = zeros(npars,3);
for nn = 1:npars
    subA = reshape(RMatrx(nn,:,:),3,3);
    rot_pf = (subA * rot')';
    rotold_pf = (subA * rotold')';
end

% Output file
fid = fopen('Fne file.txt','w');
fprintf(fid,...
    'time,delta_n,Vn_rel,Vt_rel,F1(1),F1(2),F1(3),M1(1),M1(2),M1(3)\n');

% Movie file
M = moviein(20);

% Temporal evolution starts
for tt = 1:itmax
    time = tt*dt;
    
%     fprintf('==========================\n');
    if (mod(tt,10)==0)
    fprintf('Time step tt = %d\n',tt);
    end
    
    % (1) Update patch location in the inertial frame
    
    [x_pat,y_pat,z_pat,NormVec_IF] = F_AbsPatch(xpar,ypar,zpar,x_rel,y_rel,...
    z_rel,RMatrx,NormVec,npars,NN);
    
    % (2) Calculate induced surface charge
    [sigma,E_pat]=F_InducedChrg(npars,sigma_f_scalar,...
        kappa_air,kappa_p,DeltaArea,NN,x_pat,y_pat,z_pat,NormVec_IF,...
        Ex,Ey,Ez);
    
    % (3) Calculate electrostatic force and torque
    % Note: F_pat(patch id,particle id,dimension)
    [F_par,M_par] = F_ForceTorque(sigma,E_pat,DeltaArea,NN,npars,...
    x_pat,y_pat,z_pat,xpar,ypar,zpar);
% 
%     fprintf('Fx_par = %f, Mx_par = %f\n',F_par(1,1),M_par(1,1));
%     fprintf('Fy_par = %f, My_par = %f\n',F_par(1,2),M_par(1,2));
%     fprintf('Fz_par = %f, Mz_par = %f\n',F_par(1,3),M_par(1,3));
    
    % (4) Interparticle collision
    if (npars >1)
    % Construct collision detection list
    [DetectList] = F_DetectList(npars,rpar,a,b,c);

    % Check collision pairs & identify collision point & compute contact
    % force/torques
    [CollideList,ContactPoints,vpar,rot_pf,F_par,M_par] =...
        F_Collision(npars,Elastic,Poisson,DetectList,a,b,c,...
        rpar,RMatrx,thickness,vpar,rot_pf,frict,F_par,M_par,time,fid,erest,mass);
    end % if (npars >1)
    
    % (5) Evolve ellipsoid movements
    
    % Update particle movements
    
    
    [Fold_par,vpar,vparold,rpar,rot_pf,rotold_pf,thetapar,drotxdtold,...
    drotydtold,drotzdtold] = F_Movement(npars,RMatrx,dt,mass,...
    Ixpf,Iypf,Izpf,F_par,Fold_par,vpar,vparold,rpar,...
    M_par,rot_pf,rotold_pf,thetapar,drotxdtold,drotydtold,drotzdtold);
    
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
    [x_pat,y_pat,z_pat,NormVec_IF] = F_AbsPatch(xpar,ypar,zpar,x_rel,y_rel,...
    z_rel,RMatrx,NormVec,npars,NN);
        
    figure(1);
        
    xplot = reshape(x_pat,[npars*NN,1]);
    yplot = reshape(y_pat,[npars*NN,1]);
    zplot = reshape(z_pat,[npars*NN,1]);
    scatter3(xplot,yplot,zplot,10,sigma,'filled');
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(gca,'LineWidth',1.5);
    axis equal;
    xlim([-5,5]);
    ylim([-3,3]);
    zlim([-3,3]);
    
    M(figure_id) = getframe;
    [FigA,FigMap] = rgb2ind(frame2im(getframe),256);
    if (figure_id==1)
        imwrite(FigA,FigMap,'MovieTest.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(FigA,FigMap,'MovieTest.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    figure_id = figure_id + 1;
    
    if (mod(tt,50)==0) %
    MyCheckQtn = QtnEps1(1)^2+QtnEps2(1)^2+QtnEps3(1)^2+QtnEta(1)^2;
    fprintf('Sum of quaternion square = %f\n',MyCheckQtn);
    CheckRMatrx = reshape(RMatrx(1,:,:),[3,3]);
    MyCheckRMatrx = CheckRMatrx'*CheckRMatrx;
    fprintf('[%f,%f,%f]\n',MyCheckRMatrx(1,1),MyCheckRMatrx(1,2),MyCheckRMatrx(1,3));
    fprintf('[%f,%f,%f]\n',MyCheckRMatrx(2,1),MyCheckRMatrx(2,2),MyCheckRMatrx(2,3));
    fprintf('[%f,%f,%f]\n',MyCheckRMatrx(3,1),MyCheckRMatrx(3,2),MyCheckRMatrx(3,3));
    end
    
%     % Plot the contact point
%     if (sum(CollideList)>0)
%     
%     hold on;
%     scatter3(ContactPoints(1,2,1),ContactPoints(1,2,2),...
%     ContactPoints(1,2,3),20,'r','filled');
%     scatter3(ContactPoints(2,1,1),ContactPoints(2,1,2),...
%     ContactPoints(2,1,3),20,'b','filled');
% 
%     end % if (sum(CollideList)>0)
    
end % for tt = 1:itmax 

fclose(fid);

% movie(M,1);