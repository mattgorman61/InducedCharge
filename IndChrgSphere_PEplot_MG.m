clear; close all; clc;
% INDUCED CHARGE ON A SPHERE
fprintf('INDUCED CHARGE ON A SPHERE\n\n');
% Notes for this code: https://livejohnshopkins-my.sharepoint.com/:p:/g/personal/mgorma18_jh_edu/EeEbjWCeZedBjvwOWXS5evsBl3YFRVGzGH_HL-aBKbnk0Q?e=TOGlLM
% (Slide 50)

addpath 'C:\Users\Matt Gorman\Desktop\Induced Charge MG\2021-02-11\Functions';

% Logicals
lshowNVects = true;
lshowSurfaceCharge = true;
lshowPEResults = true;

% Sphere and Medium Parameters
R = 1;
Npatches = 1500;
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



%% Electrostatic PE Plots


[U_pCharge_Norm,U_tot_Norm,netCharge] = F_getPE_Loops(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);
%[U_pCharge_Norm2,U_tot_Norm2,netCharge2] = F_getPE_Matrix(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,k_air,k_obj,epsilon_0);

%{
%Check if PE result is the same...
diffPE = U_pCharge_Norm - U_pCharge_Norm2
%}

PCx_vect = [1.05, 1.075, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5];
R_vect = PCx_vect - 1;

% kt = 0 0.025 0.1 0.4 2.5 10 40

k_air = 1; k_obj = 40;
U_kt40_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt40_Norm(i) = U_pCharge_Norm_i;
end


k_air = 1; k_obj = 10;
U_kt10_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt10_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 2.5;
U_kt2p5_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt2p5_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.4;
U_kt0p4_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p4_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.1;
U_kt0p1_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p1_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.025;
U_kt0p025_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p025_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0;
U_kt0_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0_Norm(i) = U_pCharge_Norm_i;
end


%{
% MG Results:
R_vect = [0.1, 0.2, 0.3, 0.4, 0.5];
U_vect_Norm_kt0p1 = [ 0.104330, 0.037109, 0.018965, 0.011265, 0.007305]; Net_charge_kt0p1 = [ -0.249897, -0.016342, -0.013731, -0.011961, -0.010442];
U_vect_Norm_kt0p4 = [0.056513, 0.020443, 0.010555, 0.006314, 0.004115]; Net_charge_kt0p4 = [-0.167049, -0.011802, -0.009907, -0.008604, -0.007493];
U_vect_Norm_kt2p5  = [-0.061684, -0.023304, -0.012360, -0.007533, -0.004979]; Net_charge_kt2p5 = [0.413501, 0.035193, 0.029456, 0.025395, 0.021982];
U_vect_Norm_kt10 = [-0.123970, -0.047868, -0.025782, -0.015885, -0.010582]; Net_charge_kt10 = [2.272911, 0.212290, 0.177345, 0.152331, 0.131456];
%}

% XR Result: k_med = 0.1, k_obj = 1, x_pc = 1.5R 
% Computed normalized energy is: 0.007214
% Net charge is: 0.073093





%% PE of sphere: Analytical Soln (Equation from Barros)

d_vect_anal = linspace(0.05, 0.5, 100);
k_tildaVect = [0 0.025 0.1 0.4 2.5 10 40]';
PE_vect_anal = zeros(length(d_vect_anal),length(k_tildaVect));
k_tildaVectString = strings(length(k_tildaVect),1)';

for i = 1:length(k_tildaVect)
    k_tildaVectString(i) = strcat('\kappa = ', {' '}, num2str(k_tildaVect(i)));
end

% Plot of Sphere PEs - Analytical
for j = 1:length(k_tildaVect)
    for i = 1:length(d_vect_anal)
        PE_vect_anal(i,j) = F_AnalyticalPE_ReadOnly(k_air,pcharge,k_tildaVect(j),d_vect_anal(i),R);
    end
end

if(lshowPEResults)
    
    fig3 = figure(3);
    hold on;
    
    % Note: Cell Array of Colors Requires Curly Braces!
    colorsVect = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1], [0 0 0], [1 0.2 0.1], [0.2 0.1 1]};
    
    for i = 1:length(PE_vect_anal(1,:))
        plot(d_vect_anal,PE_vect_anal(:,i),'color',colorsVect{i});
    end
    
    %{
    scatter(R_vect, U_kt0p1_Norm, 100, [1 0.7 0], '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p4_Norm, 100, [0.6 0 1], 'd', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt2p5_Norm, 100, 'g', '>', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt10_Norm, 100, 'c', 's', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt40_Norm, 100, 'r', 's', 'filled', 'MarkerEdgeColor', 'k');
    %}
    
    scatter(R_vect, U_kt0_Norm, 50, colorsVect{1}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p025_Norm, 50, colorsVect{2}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p1_Norm, 50, colorsVect{3}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p4_Norm, 50, colorsVect{4}, 's', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt2p5_Norm, 50, colorsVect{5}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt10_Norm, 50, colorsVect{6}, 's', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt40_Norm, 50, colorsVect{7}, '^', 'filled', 'MarkerEdgeColor', 'k');
    
    
    %leg3 = legend(k_tildaVectString);
    leg3 = legend([k_tildaVectString, '\kappa = 0','\kappa = 0.025','\kappa = 0.1', '\kappa = 0.4', '\kappa = 2.5', '\kappa = 10', '\kappa = 40']);
    leg3.Location = 'southeast';
    xlab3 = xlabel('Charge-surface distance (units of sphere radius)');
    xlab3.FontName = 'Times New Roman';
    ylab3 = ylabel('Normalized Potential Energy');
    ylab3.FontName = 'Times New Roman';
    set(gcf,'position',[600,100,800,700]);

end

%}

