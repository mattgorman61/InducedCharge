clear; close all; clc;
% INDUCED CHARGE ON A SPHERE
fprintf('INDUCED CHARGE ON A SPHERE\n\n');
% Notes for this code: https://livejohnshopkins-my.sharepoint.com/:p:/g/personal/mgorma18_jh_edu/EeEbjWCeZedBjvwOWXS5evsBl3YFRVGzGH_HL-aBKbnk0Q?e=TOGlLM
% (Slide 50)

% Add Function Folder to Current Path
currFolder = pwd;
% fprintf('%s',currFolder);
path1 = strcat(currFolder,'\Functions');
path2 = strcat(currFolder,'\Functions\Plots');
addpath (path);

% Logicals
lshowNVects = false;
lshowSurfaceCharge = true;
lshowPEResults = true;

% Sphere and Medium Parameters
R = 1;
Npatches = 200;
dA = 4*pi*(R^2)/Npatches;

sigma_f = zeros(Npatches,1); % Neglecting any free charges (perfect insulator?)
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



%% Discretize Spherical Surface

x = zeros(length(Npatches),1);
y = zeros(length(Npatches),1);
z = zeros(length(Npatches),1);

% Fibonacci method
gRat = (sqrt(5.0)+1.0)/2.0; % Golden Ratio
gAng = (2.0 - gRat)*(2.0*pi);
for i = 1:Npatches
    lat = asin(-1.0+2.0*double(i)/(Npatches+1));
    lon = gAng*i;
    x(i) = R*cos(lon)*cos(lat);
    y(i) = R*sin(lon)*cos(lat);
    z(i) = R*sin(lat);
end

x=x'; y=y'; z=z';

% Normal Vector:
% nVect(i,:) = nvx_i, nvy_i, nvz_i
nVect = zeros(Npatches,3); % Normal Vectors
nVect(:,1) = x/R; nVect(:,2) = y/R; nVect(:,3) = z/R;


% Plot Sphere with Normal Vectors
if(lshowNVects)
    fig1 = figure(1);
    scatter3(x,y,z,'filled','r');
    hold on;
    quiver3(x,y,z,nVect(:,1),nVect(:,2),nVect(:,3),3,'c');
    xlim([-2*R, 2*R]); ylim([-2*R, 2*R]); zlim([-2*R, 2*R]);
    axis square; 
    % tit1 = title('Sphere Normal Vectors'); tit1.FontName = 'Times New Roman';
    % tit1.FontSize = 12;
    set(gca,'LineWidth',1.5); set(gcf,'Position',[100,400,500,400]);
end

%% CALL FUNCTIONS
[sigma_b,b] = F_getSigmaB_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
[sigma_b2,b2] = F_getSigmaB_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);

%{
%Check if sigma_b result is the same...
diffSigmaB = sigma_b2 - sigma_b
diffB = b2-b
%}

sigma = sigma_b + sigma_f;

if(lshowSurfaceCharge)
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
end



%% Electrostatic PE Calculation

%{
[U_pCharge_Norm,U_tot_Norm,netCharge] = F_getPE_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);
[U_pCharge_Norm2,U_tot_Norm2,netCharge2] = F_getPE_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);


%Check if PE result is the same...
%diffPE = U_pCharge_Norm - U_pCharge_Norm2
%}

if(lshowPEResults)
    lPlotPEValid = F_Plot_PEValidation(R,x,y,z,nVect,y_pc,z_pc,pcharge,sigma_f,epsilon_0 );
end


