clear;
clc;

% =========================================================================
% Fundamental settings
filename = 'a1b1c1_2328_4652.mat'; % Geometry: Axisymmetric ellipsoid
npars = 1;%4;
dt = 0.1;%0.064; %5E-3; %1E-3
itmax = 1; % 10000
thickness = 0.05; % Preset gap in collision detection

% ========== Dimensional Paraneters in Physical Space ==========
Density_Phy = 2500; % Unit: [kg/m^3]
sigma_Phy = 0;%6.4E-6; % Unit: [C/m2]
Elastic_Phy = 1E9;%3.05E4;%1.77E8; % Unit: Pa
Ex_Phy = 4E4; % Field strength: V/m
Ey_Phy = 0; % Field strength: V/m
Ez_Phy = 0; % Field strength: V/m
Elastic_Original = 1E9; % Unit: Pa % Original elastic modulus
%
erest = 0.7; % Restitution coefficient
frict = 0.3; % Friction coefficient
Poisson = 0.33; % Poisson ratio
epsilon = 8.854E-12; % Vacuum permittivity: 8.854E-12 F/m
kappa_air = 1; % Dielectric constance of vacuum
kappa_p = 10; % Dielectric constance of particle
ReducedRatio = Elastic_Phy/Elastic_Original; % ReducedRatio of ES force

% ========== Characteristic Scales ==========
% Note: DO NOT change very often, BE CAREFUL
U_null = 1; % Velocity scale: 1 m/s
L_null = 1E-5; % Length scale: 1E-5 m = 10 um % similar to particle size
rho_null = 1000; % Density scale: 1000 kg/m^3
sigma_null = 1.6E-6; % Charge density scale: 1.6E-6 C/m2 = 10e/(um)^2
p_null = rho_null*U_null*U_null; % Pressure scale
t_null = L_null/U_null; % Time scale

% ========== Dimentionless Paraneters in Simulation ==========
% Obtained from dimensional parameters and characteristic scales

Elastic = Elastic_Phy/p_null;
rhop = Density_Phy/rho_null; % particle density
sigma_f_scalar = linspace(1,1,npars)*sigma_Phy/sigma_null;% Surface free charge density %[0,0]
EScoeff = sigma_null^2/epsilon/rho_null/U_null/U_null; % Coefficient to 
% convert dimensionless ES force to dimensionless inertia force

%==========================================================================
% Initialization
% Vectors are used to store info of all particles
% Note: vector length should equal npars

figure_id = 1; % Movie

Ex = Ex_Phy*epsilon/sigma_null; % Dimensionless field strength
Ey = Ey_Phy*epsilon/sigma_null;
Ez = Ez_Phy*epsilon/sigma_null;

% Location
xpar = 0:2.2:(npars-1)*2.2;%[2.12,-2.11,-2.11,2.11];%linspace(0,0,npars); % size=1*npars
ypar = linspace(0,0,npars);%[1.32,1.31,-1.31,-1.31];%[-1,1];
zpar = linspace(0,0,npars);%[-1,1];
rpar = [xpar;ypar;zpar]'; % size=npars*3

xrecord = zeros(itmax,npars);
yrecord = zeros(itmax,npars);
zrecord = zeros(itmax,npars);

% Velocity
vpar = zeros(npars,3); % size=npars*3
vparold = zeros(npars,3); % size=npars*3

vparxrecord = zeros(itmax,npars);
vparyrecord = zeros(itmax,npars);
vparzrecord = zeros(itmax,npars);

% Orientation angle
thetapar = zeros(npars,3);%[0,0,pi/12;0,0,-pi/12;0,0,pi/12;0,0,-pi/12];%

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

Fxrecord = zeros(itmax,npars);
Fyrecord = zeros(itmax,npars);
Fzrecord = zeros(itmax,npars);

% Quaternion rate
Eps1rateold = zeros(npars,1);
Eps2rateold = zeros(npars,1);
Eps3rateold = zeros(npars,1);
Etarateold = zeros(npars,1);

% Read patch info
[x_rel,y_rel,z_rel,DeltaArea,NN,NormVec,a,b,c]=F_GeometryEllip(filename);
mass = 4/3*rhop*pi*a*b*c; % particle mass
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
fprintf(fid,'time,overlap,vn_rel,vt_rel,id1,id2\n');

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

    if(tt == 1)
        sig_b = zeros(NN*npars,1); 
    end

    % (2) Calculate induced surface charge
    %tic;
    t1 = cputime;
    [sigma,E_pat,sig_b]=F_InducedChrg(npars,sigma_f_scalar,...
        kappa_air,kappa_p,DeltaArea,NN,x_pat,y_pat,z_pat,NormVec_IF,...
        Ex,Ey,Ez,sig_b);
    t2 = cputime - t1;
    fprintf('CPU Time = %f\n',t2);
    %tElapse = toc;
    %fprintf('Elapsed time of ES = %f\n',tElapse);
    
    % (3) Calculate electrostatic force and torque
    % Note: F_pat(patch id,particle id,dimension)
    [F_par,M_par] = F_ForceTorque(sigma,E_pat,DeltaArea,NN,npars,...
    x_pat,y_pat,z_pat,xpar,ypar,zpar);

    F_par = F_par*EScoeff;
    M_par = M_par*EScoeff; %dimensionless

%     For debug 
%     M_par = M_par*rho_null*(L_null^3)*(U_null^2); % Unit:[N*m] 
%     M_par = M_par*1E15; % Unit:[nN*um]
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
        rpar,RMatrx,thickness,vpar,rot_pf,frict,F_par,M_par,time,fid,erest,mass,ReducedRatio);
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
        
    if ((tt==1)||(mod(tt,20)==0))
%     figure(99);
    
    xplot = reshape(x_pat,[npars*NN,1]);
    yplot = reshape(y_pat,[npars*NN,1]);
    zplot = reshape(z_pat,[npars*NN,1]);
    scatter3(xplot,yplot,zplot,10,sigma,'filled');
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(gca,'LineWidth',1.5);
    axis equal;
%     xlim([-10,10]);
%     ylim([-10,10]);
%     zlim([-5,5]);
    %view([0,0,8]);
    
    M(figure_id) = getframe;
    [FigA,FigMap] = rgb2ind(frame2im(getframe),256);
    if (figure_id==1)
        imwrite(FigA,FigMap,'MovieTest.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(FigA,FigMap,'MovieTest.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    figure_id = figure_id + 1;
    end % if (mod(tt,50)==0)
    
    % For debug
%     if (mod(tt,50)==0) %
%     MyCheckQtn = QtnEps1(1)^2+QtnEps2(1)^2+QtnEps3(1)^2+QtnEta(1)^2;
%     fprintf('Sum of quaternion square = %f\n',MyCheckQtn);
%     CheckRMatrx = reshape(RMatrx(1,:,:),[3,3]);
%     MyCheckRMatrx = CheckRMatrx'*CheckRMatrx;
%     fprintf('[%f,%f,%f]\n',MyCheckRMatrx(1,1),MyCheckRMatrx(1,2),MyCheckRMatrx(1,3));
%     fprintf('[%f,%f,%f]\n',MyCheckRMatrx(2,1),MyCheckRMatrx(2,2),MyCheckRMatrx(2,3));
%     fprintf('[%f,%f,%f]\n',MyCheckRMatrx(3,1),MyCheckRMatrx(3,2),MyCheckRMatrx(3,3));
%     fprintf('x = %f, y = %f, z = %f\n',rpar(1,1),rpar(1,2),rpar(1,3));
%     fprintf('x = %f, y = %f, z = %f\n',rpar(2,1),rpar(2,2),rpar(2,3));
%     end
    
    % Record particle location
    xrecord(tt,:) = rpar(:,1)';
    yrecord(tt,:) = rpar(:,2)';
    zrecord(tt,:) = rpar(:,3)';
    
    vparxrecord(tt,:) = vpar(:,1)';
    vparyrecord(tt,:) = vpar(:,2)';
    vparzrecord(tt,:) = vpar(:,3)';
    
    Fxrecord(tt,:) = F_par(:,1)';
    Fyrecord(tt,:) = F_par(:,2)';
    Fzrecord(tt,:) = F_par(:,3)';
end % for tt = 1:itmax 

fclose(fid);
