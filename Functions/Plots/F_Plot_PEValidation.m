function [finished] = F_Plot_PEValidation(markerSize,R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma_f,epsilon_0,Ext_EField_x,Ext_EField_y,Ext_EField_z)
% PLOTS SIMULATION RESULTS AGAINST ANALYTICAL RESULTS
%   Analytical Results from Barros and Luijten 2014, Phys. Rev. Letters

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);


%% PLOT SIMULATION RESULTS
PCx_vect = [1.05, 1.075, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5];
R_vect = PCx_vect - 1;

% kt = 0 0.025 0.1 0.4 2.5 10 40

k_air = 1; k_obj = 40;
U_kt40_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt40_Norm(i) = U_pCharge_Norm_i;
end


k_air = 1; k_obj = 10;
U_kt10_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
                                     
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt10_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 2.5;
U_kt2p5_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt2p5_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.4;
U_kt0p4_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p4_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.1;
U_kt0p1_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p1_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0.025;
U_kt0p025_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt0p025_Norm(i) = U_pCharge_Norm_i;
end

k_air = 1; k_obj = 0;
U_kt0_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
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


% XR Result: k_med = 0.1, k_obj = 1, x_pc = 1.5R 
% Computed normalized energy is: 0.007214
% Net charge is: 0.073093

%}




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


    hold on;
    box on;
    
    % Note: Cell Array of Colors Requires Curly Braces!
    colorsVect = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1], [0 0 0], [1 0.2 0.1], [0.2 0.1 1]};
    
    for i = 1:length(PE_vect_anal(1,:))
        plot(d_vect_anal,PE_vect_anal(:,i),'color',colorsVect{i});
    end
    
    scatter(R_vect, U_kt0_Norm, markerSize, colorsVect{1}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p025_Norm, markerSize, colorsVect{2}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p1_Norm, markerSize, colorsVect{3}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt0p4_Norm, markerSize, colorsVect{4}, 's', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt2p5_Norm, markerSize, colorsVect{5}, '^', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt10_Norm, markerSize, colorsVect{6}, 's', 'filled', 'MarkerEdgeColor', 'k');
    scatter(R_vect, U_kt40_Norm, markerSize, colorsVect{7}, '^', 'filled', 'MarkerEdgeColor', 'k');
    
    
    %leg3 = legend(k_tildaVectString);
    %leg3 = legend([k_tildaVectString, '\kappa = 0','\kappa = 0.025','\kappa = 0.1', '\kappa = 0.4', '\kappa = 2.5', '\kappa = 10', '\kappa = 40']);
    
xlabel('\it{\bf{Charge-surface Distance (units of sphere radius)}}','fontsize',18,'fontname','Times New Roman');
%xlabel('Charge-surface Distance','fontsize',24,'fontname','Times New Roman','fontangle','Italic');
ylabel('\it{\bf{Potential Energy}}','fontsize',18,'fontname','Times New Roman');
set(gca, 'LineWidth', 2.0 );
set(gca, 'fontsize', 24.0 );
set(gca, 'XMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );
set(gca, 'YMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );
set(gca, 'XTick', 0:0.1:0.5 );
set(gca, 'YTick', -0.4:.1:0.3 );

h = findobj(gca,'Type','line');
set(h, 'Markersize', 8);
set(h, 'Linewidth', 2);

leg = legend([k_tildaVectString, '\kappa = 0','\kappa = 0.025','\kappa = 0.1', '\kappa = 0.4', '\kappa = 2.5', '\kappa = 10', '\kappa = 40']);
set(leg, 'fontsize', 16.0,'fontname','Times New Roman','fontangle','Italic','location','southeast','NumColumns',2);
%set(leg,'fontname','Times New Roman');
%set(leg,'fontangle','Italic');
axis([0 0.55 -0.4 0.3]);
set(gcf,'position',[0,0,800,600]);


finished = true;

%}
end

