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
path3 = strcat(currFolder,'\ParticleFiles');       addpath(path3);
addpath 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\Axisymmetric Particles\ParticleFiles';


% Logicals
% <<<<<<< HEAD:Validation Cases/IndChrgSphere_SingleParticle_PotentialValidation_MGUPDATED.m
lshowSpheres = false;
    lshowNVects = false;
lshowSurfaceCharge = false;
lshowPEResults = false; %  N/A for multiple spheres
lshowPotentials = true;
lshowPotentialValidation = true;
lshowForces = false;
% =======
lshowNVects = true;
lshowSurfaceCharge = true;
lshowPEResults = false; %  N/A for multiple spheres
lshowForceResults = true;
% >>>>>>> safeSpace:Axisymmetric Particles/IndChrgAxi_Mult_Master_callFunctions_20210616.m

lpCharge = false;
lEField = false;


% Sphere and Medium Parameters
R0 = 1;
% <<<<<<< HEAD:Validation Cases/IndChrgSphere_SingleParticle_PotentialValidation_MGUPDATED.m
%R = [R0,0.5*R0];
R = R0;
NpatchesSph = 1000; % Number of patches per sphere
numSpheres = 1; 
% =======
R = [R0, R0, R0];
%R = R0;
NpatchesSph = 2000; % Number of patches per sphere
numSpheres = 2; 
% >>>>>>> safeSpace:Axisymmetric Particles/IndChrgAxi_Mult_Master_callFunctions_20210616.m
Npatches = numSpheres*NpatchesSph;


% Center x,y,z coordinates for each of the spheres
dxs = [5,5,5]; dys = [0, 2.25*R0, 5*R0]; dzs = [0,0.1,0];

sigma_f = 100*ones(Npatches,1); % Neglecting any free charges (perfect insulator?)
k_obj = 1000;
k_air = .1;
k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);
%epsilon_0 = 8.85*10^(-12);
epsilon_0 = 1; % Is epsilon_0 the same as k_air?

%%{
% External E-Field NEED TO INCLUDE
Ext_EField_x = 0;
Ext_EField_y = 0;
Ext_EField_z = 0;
%}

% Point Charge Parameters
x_pcs = [1.5*R0];
y_pcs = [0];
z_pcs = [0];
%surfDists = sqrt((x_pcs-dxs).^2 + (y_pcs-dys).^2 + (z_pcs-dzs).^2) /R;
pcharge = [0];

% EVENTUALLY WANT TO MAKE CODE ABLE TO HANDLE ANY NUMBER OF POINT CHARGES



%% Discretize Spherical Surface

%{
x = zeros(Npatches,1)';
y = zeros(Npatches,1)';
z = zeros(Npatches,1)';

% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i    
nVect = zeros(Npatches,3); % Normal Vectors

% Fibonacci method
gRat = (sqrt(5.0)+1.0)/2.0; % Golden Ratio
gAng = (2.0 - gRat)*(2.0*pi);
for n = 1:numSpheres
for i = 1:NpatchesSph
    lat = asin(-1.0+2.0*double(i)/(NpatchesSph+1));
    lon = gAng*i;
    
    x(i+(n-1)*NpatchesSph) = R(n)*cos(lon)*cos(lat) + dxs(n);
    y(i+(n-1)*NpatchesSph) = R(n)*sin(lon)*cos(lat) + dys(n);
    z(i+(n-1)*NpatchesSph) = R(n)*sin(lat) + dzs(n);
    
    % Normal Vector:
    % nVect(i,:) = nvx_i, nvy_i, nvz_i
    nVect(i+(n-1)*NpatchesSph,1) = (x(i+(n-1)*NpatchesSph)-dxs(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,2) = (y(i+(n-1)*NpatchesSph)-dys(n))/R(n); 
    nVect(i+(n-1)*NpatchesSph,3) = (z(i+(n-1)*NpatchesSph)-dzs(n))/R(n);
end
end

x=x'; y=y'; z=z';
%}

%% Create Axisymmetric Spheres

% Axisymm Sphere 1
sphere1data = load('Axisymm_roughParticle4.xyz');
x1 = sphere1data(:,1); y1 = sphere1data(:,2); z1 = sphere1data(:,3);
nvx1 = sphere1data(:,1); nvy1 = sphere1data(:,2); nvz1 = sphere1data(:,3);

figure;
hold on; box on; grid on;
scatter3(x1,y1,z1,20,'filled','k');
quiver3(x1,y1,z1,nvx1,nvy1,nvz1,6,'r');
axis equal;


% <<<<<<< HEAD:Validation Cases/IndChrgSphere_SingleParticle_PotentialValidation_MGUPDATED.m
% Bad Case, Equal Angles. Hardcode so that NpatchesSph is always a perfect square
%[x,y,z,dA,dAmat,nVect,sphereID] = F_createSpheresEquiAng(R,NpatchesSph,numSpheres,dxs,dys,dzs);

% Create Ring (represents axisymmetric particle) ..... UPDATE - RING DOES
% NOT REPRESENT AXISYMMETRIC PARTICLE. FOOLISH APPROACH
%[x,y,z,dA,dAmat,nVect,sphereID] = F_createRing(R,NpatchesSph,numSpheres,dxs,dys,dzs);
% =======
% Axisymm Sphere 2
sphere2data = load('Axisymm_roughParticle2.xyz');
x2 = sphere2data(:,1); y2 = sphere2data(:,2); z2 = sphere2data(:,3);
nvx2 = sphere2data(:,1); nvy2 = sphere2data(:,2); nvz2 = sphere2data(:,3);

% figure;
% hold on; box on; grid on;
% scatter3(x2,y2,z2,20,'filled','k');
% quiver3(x2,y2,z2,nvx2,nvy2,nvz2,6,'r');
% axis equal;


% Axisymm Sphere 3
sphere3data = load('Axisymm_roughParticle3.xyz');
x3 = sphere3data(:,1); y3 = sphere3data(:,2); z3 = sphere3data(:,3);
nvx3 = sphere3data(:,1); nvy3 = sphere3data(:,2); nvz3 = sphere3data(:,3);

% figure;
% hold on; box on; grid on;
% scatter3(x3,y3,z3,20,'filled','k');
% quiver3(x3,y3,z3,nvx3,nvy3,nvz3,6,'r');
% axis equal;


x = [x1; x1 + 0];
nvx = [nvx1;nvx1];
y = [y1; y1 + 2200];
nvy = [nvy1;nvy1];
z = [z1; z1 + 0];
nvz = [nvz1;nvz1];

nVect = [nvx,nvy,nvz];
% >>>>>>> safeSpace:Axisymmetric Particles/IndChrgAxi_Mult_Master_callFunctions_20210616.m


% Plot Sphere with (or without) Normal Vectors
if(lshowSpheres)
    figure();
% <<<<<<< HEAD:Validation Cases/IndChrgSphere_SingleParticle_PotentialValidation_MGUPDATED.m
    F_Plot_NormVectors(R,x,y,z,nVect,lshowNVects);
% =======
%     axis equal;
    F_Plot_NormVectors(R0,x,y,z,nVect,lshowNVects);
% >>>>>>> safeSpace:Axisymmetric Particles/IndChrgAxi_Mult_Master_callFunctions_20210616.m
end


%% CALL FUNCTIONS
%Multiple Spheres:
[sigma_b,b] = F_getSigmaB_Mult_Matrix(dAmat,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);

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
    F_Plot_sigmaB(x,y,z,x_pcs,y_pcs,z_pcs,sigma_b,lpCharge,lEField,Ext_EField_x,Ext_EField_y,Ext_EField_z);
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

%sa;lkfja;slkdfsj



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

%% Electrostatic Potential Calculation

if(lshowPotentials)
    figure();
    F_Plot_Potential(R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0);
end

if(lshowPotentialValidation)
    figure();
    F_Plot_PotentialValidation(R,x,y,z,nVect,x_pcs,y_pcs,z_pcs,pcharge,sigma,epsilon_0,Ext_EField_x,Ext_EField_y,Ext_EField_z);
end




%% Electrostatic Force Calculation

if(lshowForces)
    figure();
% <<<<<<< HEAD:Validation Cases/IndChrgSphere_SingleParticle_PotentialValidation_MGUPDATED.m
    F_Plot_ForceVectors(R,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0);
end


%}
% =======
    F_Plot_ForceVectors(R0,x,y,z,dA,nVect,x_pcs,y_pcs,z_pcs,pcharge,lpCharge,lEField,sigma_b,sigma_f,k_air,k_obj,epsilon_0);
% end


% >>>>>>> safeSpace:Axisymmetric Particles/IndChrgAxi_Mult_Master_callFunctions_20210616.m
