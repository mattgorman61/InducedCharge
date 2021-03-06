clear; close all; clc;
% INDUCED CHARGE ON A SPHERE
fprintf('INDUCED CHARGE ON A SPHERE\n\n');
% Notes for this code: https://livejohnshopkins-my.sharepoint.com/:p:/g/personal/mgorma18_jh_edu/EeEbjWCeZedBjvwOWXS5evsBl3YFRVGzGH_HL-aBKbnk0Q?e=TOGlLM
% (Slide 50)

% Add Function Folder to Current Path
currFolder = pwd;
% fprintf('%s',currFolder);
path1 = strcat(currFolder,'\..\Functions');        addpath(path1);
path2 = strcat(currFolder,'\..\Functions\Plots');  addpath(path2);


% Logicals
lshowNVects = false;
lshowSurfaceCharge = true;
lshowPEResults = true;
lshowForceResults = true;


% Sphere and Medium Parameters
R = 1;
NpatchesSph = 200;
dA = 4*pi*(R^2)/NpatchesSph;
numSpheres = 1;

sigma_f = zeros(NpatchesSph,1); % Neglecting any free charges (perfect insulator?)
k_obj = 1;
k_air = .1;
k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);
%epsilon_0 = 8.85*10^(-12);
epsilon_0 = 1;

% Point Charge Parameters
x_pc = 1.5*R;
surfDist = (x_pc-R)/R;
y_pc = 0;
z_pc = 0;
pcharge = -1;


% External Field Parameters
Ext_EField_x = 0;
Ext_EField_y = 0;
Ext_EField_z = 0;


%% Discretize Spherical Surface

x = zeros(NpatchesSph,1);
y = zeros(NpatchesSph,1);
z = zeros(NpatchesSph,1);

% Fibonacci method
gRat = (sqrt(5.0)+1.0)/2.0; % Golden Ratio
gAng = (2.0 - gRat)*(2.0*pi);
for i = 1:NpatchesSph
    lat = asin(-1.0+2.0*double(i)/(NpatchesSph+1));
    lon = gAng*i;
    x(i) = R*cos(lon)*cos(lat);
    y(i) = R*sin(lon)*cos(lat);
    z(i) = R*sin(lat);
end

x=x'; y=y'; z=z';

% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i
nVect = zeros(NpatchesSph,3); % Normal Vectors
nVect(:,1) = x/R; nVect(:,2) = y/R; nVect(:,3) = z/R;


% Plot Sphere with Normal Vectors
if(lshowNVects)
    figure();
    F_Plot_NormVectors(R,x,y,z,nVect);
end

%% CALL FUNCTIONS
% [sigma_b,b] = F_getSigmaB_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
[sigma_b2,b2] = F_getSigmaB_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
%R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z
%{
%Check if sigma_b result is the same...
diffSigmaB = sigma_b2 - sigma_b
diffB = b2-b
%}

sigma = sigma_b2 + sigma_f;

if(lshowSurfaceCharge)
    figure();
    F_Plot_sigmaB(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);
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
    lPlotPEValid = F_Plot_PEValidation(R,x,y,z,nVect,y_pc,z_pc,pcharge,sigma_f,epsilon_0 );
end

%% Electrostatic Force Calculation

if(lshowForceResults)
    figure();
    F_Plot_ForceVectors(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_b,sigma_f,k_air,k_obj,epsilon_0);
end


