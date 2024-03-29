clear; close all; clc;
% INDUCED CHARGE ON A SPHERE
fprintf('INDUCED CHARGE ON AXISYMMETRIC SPHERES\n\n');
% Notes for this code: https://livejohnshopkins-my.sharepoint.com/:p:/g/personal/mgorma18_jh_edu/EeEbjWCeZedBjvwOWXS5evsBl3YFRVGzGH_HL-aBKbnk0Q?e=TOGlLM
% (Slide 50)

% Add Function Folder to Current Path
currFolder = pwd;
% fprintf('%s',currFolder);
path1 = strcat(currFolder,'\..\Functions');        addpath(path1);
path2 = strcat(currFolder,'\..\Functions\Plots');  addpath(path2);
path3 = strcat(currFolder,'\..\Rough Particles');  addpath(path3);

% Logicals
lshowNVects = true;
lshowSurfaceCharge = true;
lshowPEResults = false; %  N/A for multiple spheres
lshowForceResults = false;

lpCharge = false;
lEField = false;


% Sphere and Medium Parameters
R0 = 1;
R = [R0, R0, R0];
%R = R0;
NpatchesSph = 3934; % Number of patches per sphere. Have to change every time distmesh settings are changed.
numSpheres = 2; 
Npatches = numSpheres*NpatchesSph;


% x,y,z coordinates for centers of each of the spheres
dxs = [0,0,0]; dys = [0, 2.3*R0, 4.5*R0]; dzs = [0,0,0]; % For Tip-to-Tip configuration
% dxs = [0,0,0]; dys = [0, 2.25*R0, 4.5*R0]; dzs = [0, 0.25*R0, 0.5*R0]; % For Mountain-Valley configuration

sigma_f = [100000*ones(Npatches/2,1); 100000*ones(Npatches/2,1)]; % Neglecting any free charges (perfect insulator?)
k_obj = 1;
k_air = .1;
k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);
%epsilon_0 = 8.85*10^(-12);
epsilon_0 = 1;

%%{
% External E-Field NEED TO INCLUDE
Ext_EField_x = 0;
Ext_EField_y = 0;
Ext_EField_z = 0;
%}

% Point Charge Parameters
x_pcs = [15*R0];
y_pcs = [0];
z_pcs = [0];
%surfDists = sqrt((x_pcs-dxs).^2 + (y_pcs-dys).^2 + (z_pcs-dzs).^2) /R;
pcharge = [0];

% EVENTUALLY WANT TO MAKE CODE ABLE TO HANDLE ANY NUMBER OF POINT CHARGES


%% Create Particles: Spheres or Axisymmetric Spheres
x = zeros(Npatches,1);
y = zeros(Npatches,1);
z = zeros(Npatches,1);
dA = zeros(Npatches,1);
nVect = zeros(Npatches,3);
nvMags = zeros(Npatches,1);

for n = 1:numSpheres

    A = 0.25; alpha = 12; elemSize = 0.1;
    
%     [x_n,y_n,z_n,dA_n,nVect_n,sphereID_n] = F_createEllipPar_distmesh(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);
    [x_n,y_n,z_n,dA_n,V_n,nVect_n] = F_createAxiSymmPar_distmesh(R(n),dxs(n),dys(n),dzs(n),A,alpha,elemSize);
%     [x_n,y_n,z_n,dA_n,dAmat,nVect_n,sphereID_n] = F_createSingSphere(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);
%     [x_n,y_n,z_n,dA_n,nVect_n,sphereID_n] = F_createAxisymmSphere(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);
%     [x_n,y_n,z_n,dA_n,nVect_n,nvMags_n,sphereID_n] = F_createAxisymmSphere_ThetaSymm(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);
%     [x_n,y_n,z_n,dA_n,nVect_n,sphereID_n] = F_createAxisymmSphere_PhiSymm(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);
%     [x_n,y_n,z_n,dA_n,nVect_n,sphereID_n] = F_createSimpleSphere(R(n),NpatchesSph,dxs(n),dys(n),dzs(n),n);

    x((n-1)*NpatchesSph + 1: n*NpatchesSph) = x_n;
    y((n-1)*NpatchesSph + 1: n*NpatchesSph) = y_n;
    z((n-1)*NpatchesSph + 1: n*NpatchesSph) = z_n;
    dA((n-1)*NpatchesSph + 1: n*NpatchesSph) = dA_n;
    nVect((n-1)*NpatchesSph + 1: n*NpatchesSph,1) = nVect_n(:,1);
    nVect((n-1)*NpatchesSph + 1: n*NpatchesSph,2) = nVect_n(:,2);
    nVect((n-1)*NpatchesSph + 1: n*NpatchesSph,3) = nVect_n(:,3);
%     nvMags((n-1)*NpatchesSph + 1: n*NpatchesSph) = nvMags_n;
    
end

dAmat = repmat(dA',Npatches,1);

% Plot Sphere with Normal Vectors
if(lshowNVects)
    figure();
%     axis equal;
    F_Plot_NormVectors(R0,x,y,z,nVect,lshowNVects);
end

p = [x,y,z];
nvMags = sqrt(nVect(:,1).*nVect(:,1) + nVect(:,2).*nVect(:,2) + nVect(:,3).*nVect(:,3));
% save('AxisymmPar_PhiSymm_300_E','dA','a','b','c','p','nVect' );
% save('ParTEST_distmesh_576','p','nVect','dA');

% xc,yc,zc,DeltaArea,N,NormVec,a,b,c

%% CALL FUNCTIONS
%Multiple Spheres:
[sigma_b,b] = F_getSigmaB_Mult_Matrix(numSpheres,NpatchesSph,dAmat,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);

%Single Sphere:
%[sigma_b,b] = F_getSigmaB_Matrix(R,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);

%{
%Check if sigma_b result is the same...
diffSigmaB = sigma_b2 - sigma_b
diffB = b2-b
%}

sigma = sigma_b + sigma_f;

if(lshowSurfaceCharge)
    figure();
    patchSize = 50;
    F_Plot_sigmaB(patchSize,x,y,z,x_pcs,y_pcs,z_pcs,sigma_b,lpCharge,lEField,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    
    %{
    fig2 = figure(2);
    scatter3(x,y,z,12,sigma_b,'filled');
    hold on; 
    scatter3(x_pc,y_pc,z_pc,12,'filled','k');
    xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);
    colorbar;
    axis square;    
    tit2 = title('Surface Bound Charge'); tit2.FontSize = 12; 
    tit2.FontName = 'Times New Roman';
    set(gca,'LineWidth',1.5); set(gcf,'Position',[100,100,500,400]);
    view(35,20);
    %}
end



%% Electrostatic PE Calculation

%{
[U_pCharge_Norm,U_tot_Norm,netCharge] = F_getPE_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);
[U_pCharge_Norm2,U_tot_Norm2,netCharge2] = F_getPE_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);


%Check if PE result is the same...
%diffPE = U_pCharge_Norm - U_pCharge_Norm2
%}

if(lshowPEResults)
    figure();
    F_Plot_PEValidation(R,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,epsilon_0,Ext_EField_x,Ext_EField_y,Ext_EField_z );
end

%% Electrostatic Force Calculation

if(lshowForceResults)
    figure();
    F_Plot_ForceVectors(numSpheres,NpatchesSph,R0,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,dxs,dys,dzs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0);
end


