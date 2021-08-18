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

fileList1 = dir(fullfile(pwd,'*.mat'));
numPars = length(fileList1);
NpatchesPar = zeros(numPars,1);
for i = 1:length(fileList1)
    pardata_i = load(fileList1(i).name);
    ParticleData.x{i} = pardata_i.p(:,1);
    ParticleData.y{i} = pardata_i.p(:,2);
    ParticleData.z{i} = pardata_i.p(:,3);
    ParticleData.dA{i} = pardata_i.dA;
    ParticleData.V{i} = pardata_i.V;
    ParticleData.nvect{i} = pardata_i.nVect;
    
    NpatchesPar(i) = length(ParticleData.x{1,i});
end

Npatches = sum(NpatchesPar);


% x,y,z coordinates for centers of each of the particles
R0 = 1;
dxs = [0,0,10*R0]; dys = [0, 2.2*R0, 10*R0]; dzs = [0,-10*R0,-25*R0]; % For Tip-to-Tip configuration
% dxs = [0,0,0]; dys = [0, 2.25*R0, 4.5*R0]; dzs = [0, 0.25*R0, 0.5*R0]; % For Mountain-Valley configuration

sigma_f = [1*ones(NpatchesPar(1),1); -1*ones(NpatchesPar(2),1)]; % Neglecting any free charges (perfect insulator?)
k_obj = 1;
k_air = .1;
k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);
%epsilon_0 = 8.85*10^(-12);
epsilon_0 = 1;

%%{
% External E-Field 
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


%% Create Particles: Spheres or Axisymmetric Spheres
x = zeros(Npatches,1);
y = zeros(Npatches,1);
z = zeros(Npatches,1);
dA = zeros(Npatches,1);
nVect = zeros(Npatches,3);
nvMags = zeros(Npatches,1);

for n = 1:numPars
  
    if (n == 1)
        x(1:NpatchesPar(1)) = ParticleData.x{n} + dxs(n);
        y(1:NpatchesPar(1)) = ParticleData.y{n} + dys(n);
        z(1:NpatchesPar(1)) = ParticleData.z{n} + dzs(n);
        
        dA(1:NpatchesPar(1)) = ParticleData.dA{n};
        nVect(1:NpatchesPar(1),:) = ParticleData.nvect{1,n};
    else
        index = sum(NpatchesPar(1:n-1));
        x(index + 1 : index+NpatchesPar(n)) = ParticleData.x{1,n} + dxs(n);
        y(index + 1 : index+NpatchesPar(n)) = ParticleData.y{1,n} + dys(n);
        z(index + 1 : index+NpatchesPar(n)) = ParticleData.z{1,n} + dzs(n);
        dA(index + 1 : index+NpatchesPar(n)) = ParticleData.dA{n};
        nVect(index + 1 : index+NpatchesPar(n),:) = ParticleData.nvect{1,n};
    end
    
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
% Multiple Particles:
[sigma_b,b] = F_getSigmaB_Arb_Matrix(numPars,NpatchesPar,dAmat,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);

% Single Sphere:
% [sigma_b,b] = F_getSigmaB_Matrix(R,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);

sigma = sigma_b + sigma_f;

if(lshowSurfaceCharge)
    figure();
    patchSize = 30;
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


