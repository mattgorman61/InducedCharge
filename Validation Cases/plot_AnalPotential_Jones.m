clc; clear all; close all;

R = 1;
pcharge = -1;
x_pc = 1.5*R;
y_pc = 0;
z_pc = 0;
k_air = 1;
k_sphere = [0.1, 0.5, 1, 2.5, 10, 50];
N = 20;


plotVect = zeros(length(k_sphere));
legstring = strings(length(k_sphere),1);

figure();
hold on;
box on;
    
for i = 1:length(k_sphere)    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,x_pc,y_pc,z_pc,k_air,k_sphere(i),N);
    p_i = plot(theta_i,spherePot_i);
    plotVect(i) = p_i;
    legstring(i) = strcat('\kappa_{s} / \kappa_{air} = ', num2str(k_sphere(i)/k_air,'%0.2f') );
        
end
legend(legstring);
title('N = ', N);

figure();
hold on;
box on;
Nvect = [1 5 10 20 50 55 ];
legstring2 = strings(length(Nvect),1);

for i = 1:length(Nvect)    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,x_pc,y_pc,z_pc,k_air,k_sphere(4),Nvect(i));
    p_i = plot(theta_i,spherePot_i);
    plotVect(i) = p_i;
    legstring2(i) = strcat('N = ', num2str(Nvect(i)));
end
title_str = ['\kappa_s / \kappa_{air} = ', num2str(k_sphere(4)/k_air),' ', 'd/R = ', num2str(x_pc/R - 1)];
title (title_str);
legend(legstring2);

figure();
hold on;
box on;
x_pcVect = R * [1.05, 1.1, 1.2, 1.3, 1.4, 1.5];
N=20;
    
for i = 1:length(x_pcVect)    
    [theta_i,spherePot_i] = F_get_AnalPotentialSphere_Jones(R,pcharge,x_pcVect(i),y_pc,z_pc,k_air,k_sphere(i),N);
    p_i = plot(theta_i,spherePot_i);
    plotVect(i) = p_i;
    legstring(i) = strcat('d/R = ', num2str(x_pcVect(i)/R - 1,'%0.2f') );
        
end
legend(legstring);
%title('N = ', N);








