function [finished] = F_Plot_PotentialValidation(R,x,y,z,nVect,x_pc,y_pc,z_pc,pcharge,sigma,epsilon_0,Ext_EField_x,Ext_EField_y,Ext_EField_z)
% PLOTS SIMULATION RESULTS AGAINST ANALYTICAL RESULTS
%   Analytical Results are from Jones, Electromechanics of Particles

% Add Functions Folder to the path
currDir = pwd;
%fprintf('%s',currDir);
idcs   = strfind(currDir,'\');
newdir = currDir(1:idcs(end)-1);
addpath (newdir);


%% PLOT SIMULATION RESULTS
PCx_vect = [1.05, 1.1, 1.2, 1.3, 1.4, 1.5];
R_vect = PCx_vect - 1;
k_air = 1;
k_obj = [0.1, 0.5, 1, 2.5, 10, 50];
N = 20;
dA = zeros(length(x));

% Note: Cell Array of Colors Requires Curly Braces!
colorsVect = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [0 1 1], [1 0 1], [0 0 0], [1 0.2 0.1], [0.2 0.1 1]};


hold on;
legstring2 = strings(length(PCx_vect),1);

z0indices = find(abs(z)<0.05)
x0 = x(z0indices);
y0 = y(z0indices);
z0 = z(z0indices);

%{ 
%Original: Ring Geometry
for i = 1:length(PCx_vect)
    ii = mod(i,length(colorsVect));
    colorsVect{ii};
        
    [phi,phi_0,phi_norm,theta] = F_getPotentials_Matrix(R,x,y,z,dA,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma,k_air,k_obj(4),epsilon_0);
    scatter(theta(1:uint8(length(theta)/8)), phi_norm(1:uint8(length(phi_norm)/8)), 50, colorsVect{ii} );
    legstring2(i) = ['d/R = ', ' ', num2str(R_vect(i))];
end
%}

%%{
% Updated: extract z0 indices from sphere geometry
for i = 1:length(PCx_vect)
    ii = mod(i,length(colorsVect));
    colorsVect{ii};
        
    [phi,phi_0,phi_norm,theta] = F_getPotentials_Matrix(R,x,y,z,dA,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma,k_air,k_obj(4),epsilon_0);
    scatter(theta, phi_norm, 50, colorsVect{ii} );
    legstring2(i) = ['d/R = ', ' ', num2str(R_vect(i))];
end
%}

%{
k_air = 1; k_obj = 40;
U_kt40_Norm = zeros(length(PCx_vect),1);
for i = 1:length(PCx_vect)
    [sigma_b_i] = F_getSigmaB_Matrix(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_f,k_air,k_obj,Ext_EField_x,Ext_EField_y,Ext_EField_z);
    sigma_i = sigma_b_i + sigma_f;
    [U_pCharge_Norm_i] = F_getPE_Loops(R,x,y,z,nVect,PCx_vect(i),y_pc,z_pc,pcharge,sigma_i,k_air,k_obj,epsilon_0);
    U_kt40_Norm(i) = U_pCharge_Norm_i;
end

%}






%% PE of sphere: Analytical Soln (Equation from Jones)

plotVect = zeros(length(k_obj));
legstring = strings(length(k_obj),1);

%figure();
hold on;
box on;
N=20;
    
for i = 1:length(PCx_vect)  
    ii = mod(i,length(colorsVect));
    colorsVect{ii};
    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,PCx_vect(i),y_pc,z_pc,k_air,k_obj(4),N);
    
    p_i = plot(theta_i,spherePot_i,'Color',colorsVect{ii});
    plotVect(i) = p_i;
    legstring(i) = ['d/R = ', num2str(R_vect(i)) ];
        
end

legstringAll = [legstring2,legstring];

legend(legstringAll,'location','southeast','NumColumns',2);
%title('N = ', N); 
    
xlabel('\it{\bf{Charge-surface Distance (units of particle radius)}}','fontsize',24,'fontname','Times New Roman');
%xlabel('Charge-surface Distance','fontsize',24,'fontname','Times New Roman','fontangle','Italic');
ylabel('\it{\bf{Potential (\psi)}}','fontsize',24,'fontname','Times New Roman');
set(gca, 'LineWidth', 2.0 );
set(gca, 'fontsize', 16.0 );
set(gca, 'XMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );
set(gca, 'YMinorTick', 'on');
set(gca, 'Ticklength', [0.02;0.01] );
%set(gca, 'XTick', 0:0.1:0.5 );
%set(gca, 'YTick', -0.4:.1:0.3 );

h = findobj(gca,'Type','line');
set(h, 'Markersize', 8);
set(h, 'Linewidth', 2);

%leg = legend([k_tildaVectString, '\kappa = 0','\kappa = 0.025','\kappa = 0.1', '\kappa = 0.4', '\kappa = 2.5', '\kappa = 10', '\kappa = 40']);
%set(leg, 'fontsize', 16.0,'fontname','Times New Roman','fontangle','Italic','location','southeast','NumColumns',2);
%set(leg,'fontname','Times New Roman');
%set(leg,'fontangle','Italic');
%axis([0 0.55 -0.4 0.3]);
set(gcf,'position',[0,0,800,600]);
xlim([0 1]);

% Display influence of k_obj/k_air on solution results
%{
figure();
hold on;
box on;
    
for i = 1:length(k_obj)    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,x_pc,y_pc,z_pc,k_air,k_obj(i),N);
    p_i = plot(theta_i,spherePot_i);
    plotVect(i) = p_i;
    legstring(i) = strcat('\kappa_{s} / \kappa_{air} = ', num2str(k_obj(i)/k_air,'%0.2f') );
        
end

legend(legstring);
title('N = %g', N);
%}

% Display influence of number of terms in the partial sum on solution results
%{
figure();
hold on;
box on;
Nvect = [1 5 10 20 50 55 ];
legstring2 = strings(length(Nvect),1);

for i = 1:length(Nvect)    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,x_pc,y_pc,z_pc,k_air,k_obj(4),Nvect(i));
    p_i = plot(theta_i,spherePot_i);
    plotVect(i) = p_i;
    legstring2(i) = strcat('N = ', num2str(Nvect(i)));
end
title_str = ['\kappa_s / \kappa_{air} = ', num2str(k_obj(4)/k_air),' ', 'd/R = ', num2str(x_pc/R - 1)];
title (title_str);
legend(legstring2);
%}


finished = true;

%}

end

