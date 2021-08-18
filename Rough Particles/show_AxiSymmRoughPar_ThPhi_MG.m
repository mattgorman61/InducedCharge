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


% Logicals
lshowNVects = false;
lshowSurfaceCharge = true;
lshowPEResults = false; %  N/A for multiple spheres
lshowForceResults = true;

lpCharge = false;
lEField = false;

% Medium Parameters
k_obj = 1;
k_air = .1;
k_tilda = k_obj/k_air; k_delta = k_air - k_obj; k_bar = 0.5*(k_air + k_obj);
%epsilon_0 = 8.85*10^(-12);
epsilon_0 = 1;



%% Create Particles: Spheres or Axisymmetric Spheres

% Sphere Parameters
R0 = 1;
R = [R0];
numSpheres = 1; 
NpatchesSph = 5000;
dxs = [0]; dys = [0]; dzs = [0];  % x,y,z coordinates for sphere center
[x_n,y_n,z_n,dA_n,V_n,nVect_n] = F_createAxisymmSphere(R,NpatchesSph,dxs,dys,dzs,1);

NpatchesSph = length(x_n);
Npatches = numSpheres*NpatchesSph;

x = x_n;
y = y_n;
z = z_n;
dA = dA_n;
V = V_n;
nVect(:,1) = nVect_n(:,1);
nVect(:,2) = nVect_n(:,2);
nVect(:,3) = nVect_n(:,3);

sigma_f = [1000*ones(Npatches/2,1); 1000*ones(Npatches/2,1)]; % Neglecting any free charges (perfect insulator?)


%     nvMags((n-1)*NpatchesSph + 1: n*NpatchesSph) = nvMags_n;


dAmat = repmat(dA',Npatches,1);

% Plot Sphere with Normal Vectors

figure();
F_Plot_NormVectors(R0,x,y,z,nVect,lshowNVects);


p = [x,y,z];
nvMags = sqrt(nVect(:,1).*nVect(:,1) + nVect(:,2).*nVect(:,2) + nVect(:,3).*nVect(:,3));
% save('AxisymmPar_PhiSymm_300_E','dA','a','b','c','p','nVect' );

filenameStr = ['Par_distmesh_',num2str(length(x))];
uinp = input(['Would you like to save this particle as ', filenameStr ,'.mat (y/n)?\n   '],'s');
if(strcmp(uinp,'y') == 1)
    save(filenameStr,'p','nVect','dA','V');
end

% xc,yc,zc,DeltaArea,N,NormVec,a,b,c

