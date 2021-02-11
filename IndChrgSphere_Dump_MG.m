clear; close all; clc;
% INDUCED CHARGE ON A SPHERE
fprintf('INDUCED CHARGE ON A SPHERE\n\n');
% Notes for this code: https://livejohnshopkins-my.sharepoint.com/:p:/g/personal/mgorma18_jh_edu/EeEbjWCeZedBjvwOWXS5evsBl3YFRVGzGH_HL-aBKbnk0Q?e=TOGlLM
% (Slide 50)
% THIS FILE DOES NOT CALL ANY FUNCTIONS. ALL THE CODE IS DUMPED HERE.

currFolder = pwd;
% fprintf('%s',currFolder);
path = strcat(currFolder,'\Functions');
addpath (path);


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
k_tilda = k_obj/k_air; k_delta = k_air - k_obj;
k_bar = 0.5*(k_air + k_obj);
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

% Normal Vector Matrix:
% nVectM(:,:,1) = nVX1, nVX2, nVX3, ... (repeated for each row)
% nVectM(:,;,2) = nVY1, nVY2, nVY3, ... (repeated for each row)
% etc.
nVectM = zeros(Npatches, Npatches,3);
for i = 1:3
nVectM(:,:,i) = repmat(nVect(:,i),1,Npatches);
end


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

%% FOR-LOOP APPROACH
%%{
A1 = zeros(Npatches);
B1 = zeros(Npatches);

% Patch contributions
for i = 1:Npatches
for j = 1:Npatches
    
    dx = x(i)-x(j); dy = y(i)-y(j); dz = z(i)-z(j);
    dr = sqrt(dx^2 + dy^2 + dz^2);
    
    if i==j
        A1(i,j) = k_bar;
        B1(i,j) = 1-k_bar;
    else
        A1(i,j) = k_delta*dA/4/pi*(dx*nVect(i,1) + dy*nVect(i,2) + dz*nVect(i,3))/(dr^3);
        B1(i,j) = -k_delta*dA/4/pi*(dx*nVect(i,1) + dy*nVect(i,2) + dz*nVect(i,3))/(dr^3);
    end
end
end

% Point charge contribution
b_pc = zeros(Npatches,1);
for i = 1:Npatches
    dx_pc = x(i)-x_pc; dy_pc = y(i)-y_pc; dz_pc = z(i)-z_pc;
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    
    b_pc(i) = (dx_pc*nVect(i,1) + dy_pc*nVect(i,2) + dz_pc*nVect(i,3))/(dr_pc^3);
end

b1 = B1*sigma_f - k_delta*pcharge/4/pi*(b_pc);
sigma_b = A1\b1;
%}

%% MATRIX APPROACH

%%{

% Patch-to-patch location differences matrix
% ppld(:,:,1) = x1-x1, x1-x2, x1-x3, ... 
%               x2-x1, x2-x2, x2-x3, ... xi-xj
ppld = zeros(Npatches,Npatches,3); 
ppld(:,:,1) = x - x'; ppld(:,:,2) = y - y'; ppld(:,:,3) = z - z';
rpp = sqrt(ppld(:,:,1).*ppld(:,:,1) + ppld(:,:,2).*ppld(:,:,2) + ...
    ppld(:,:,3).*ppld(:,:,3));

% Dot product term of A
M = (ppld(:,:,1).*nVectM(:,:,1) + ppld(:,:,2).*nVectM(:,:,2) + ...
    ppld(:,:,3).*nVectM(:,:,3))./(rpp.^3);
M(isnan(M) | isinf(M)) = 0;

A2 = k_bar*eye(Npatches) + dA*k_delta/4/pi*M;

% pCharge-to-patch location differences vector
% pcpld(:,:) = xpc-x1, ypc-y1, zpc-z1
%              xpc-x2, ypc-y2, zpc-z2
%              etc.
pcpld = zeros(Npatches,3);
pcpld(:,1) = x-x_pc; pcpld(:,2) = y-y_pc; pcpld(:,3) = z-z_pc;
rpcp = sqrt(pcpld(:,1).*pcpld(:,1) + pcpld(:,2).*pcpld(:,2) + ...
    pcpld(:,3).*pcpld(:,3));

% nVect 
% nvx(i,:) = nvx_i, nvy_i, nvz_i (Repeat for each row i)

% Pcharge-to-patch dot product matrix
Mpc = diag(pcpld*nVect')./(rpcp.^3);
%Mpc = 
Mpc(isnan(Mpc) | isinf(Mpc)) = 0;

b2 = ((1-k_bar)*eye(Npatches) - k_delta*dA/4/pi*M)*sigma_f - k_delta/4/pi*pcharge*Mpc;
sigma_b2 = A2\b2;

% Compare Results: Matrix vs. Loops
%{
diffSigmaB = sigma_b2 - sigma_b
diffB = b2-b1
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
        PE_vect_anal(i,j) = F_getAnalyticalPEResults_Barros_ReadOnly(k_air,pcharge,k_tildaVect(j),d_vect_anal(i),R);
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

% OLD CODE, BEFORE 2/11/2021
%{

U_tot_patches = 0;
U_tot_pcharge = 0;

for i = 1:Npatches
% Patch-patch contribution
for j = 1:Npatches
    if (i~=j)
    dx = x(i)-x(j); dy = y(i)-y(j); dz = z(i)-z(j);
    dr = sqrt(dx^2 + dy^2 + dz^2);
    U_tot_patches = U_tot_patches + (sigma(i)*(dA^2))/(dr);
    end
end
% Pcharge-patch contribution
    dx_pc = x(i)-x_pc;  dy_pc = y(i)-y_pc; dz_pc = z(i)-z_pc;
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    U_tot_patches = U_tot_patches + sigma(i)*dA*pcharge/(dr_pc);
end

% Patch-pcharge contribution
for i = 1:Npatches
    dx_pc = x_pc-x(i);  dy_pc = y_pc-y(i); dz_pc = z_pc-z(i);
    dr_pc = sqrt(dx_pc^2 + dy_pc^2 + dz_pc^2);
    U_tot_pcharge = U_tot_pcharge + sigma(i)*dA*pcharge/(dr_pc);
end

U_tot = 0.5*(1/4/pi/epsilon_0)*(U_tot_patches + U_tot_pcharge);
U0 = (pcharge^2)/epsilon_0/k_air/R;
%U_tot_Norm = U_tot/U0;
%netCharge = sum(sigma);

U_tot_Norm =  (0.5/4/pi/epsilon_0) * U_tot_pcharge / U0;
netCharge = sum(sigma);

fprintf('RESULTS:\n')
fprintf('k_tilda: %.5f\n', k_tilda);
fprintf('Surface Distance: %f\n', surfDist);
fprintf('Normalized PE: %f\n', U_tot_Norm);
fprintf('Net Charge: %f\n\n\n', netCharge);

% MG Results:
R_vect = [0.1, 0.2, 0.3, 0.4, 0.5];
U_vect_Norm_kt0p1 = [ 0.104330, 0.037109, 0.018965, 0.011265, 0.007305]; Net_charge_kt0p1 = [ -0.249897, -0.016342, -0.013731, -0.011961, -0.010442];
U_vect_Norm_kt0p4 = [0.056513, 0.020443, 0.010555, 0.006314, 0.004115]; Net_charge_kt0p4 = [-0.167049, -0.011802, -0.009907, -0.008604, -0.007493];
U_vect_Norm_kt2p5  = [-0.061684, -0.023304, -0.012360, -0.007533, -0.004979]; Net_charge_kt2p5 = [0.413501, 0.035193, 0.029456, 0.025395, 0.021982];
U_vect_Norm_kt10 = [-0.123970, -0.047868, -0.025782, -0.015885, -0.010582]; Net_charge_kt10 = [2.272911, 0.212290, 0.177345, 0.152331, 0.131456];

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
    for i = 1:length(PE_vect_anal(1,:))
        plot(d_vect_anal,PE_vect_anal(:,i));
    end
    
    scatter(R_vect, U_vect_Norm_kt0p1, 100, [1 0.7 0], '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_vect_Norm_kt0p4, 100, [0.6 0 1], 'd', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_vect_Norm_kt2p5, 100, 'g', '>', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_vect_Norm_kt10, 100, 'c', 's', 'filled', 'MarkerEdgeColor', 'k');
    
    %leg3 = legend(k_tildaVectString);
    leg3 = legend([k_tildaVectString, '\kappa = 0.1', '\kappa = 0.4', '\kappa = 2.5', '\kappa = 10']);
    leg3.Location = 'southeast';
    xlab3 = xlabel('Charge-surface distance (units of sphere radius)');
    ylab3 = ylabel('Normalized Potential Energy');
    set(gcf,'position',[600,100,800,600]);

end

%}
%}
