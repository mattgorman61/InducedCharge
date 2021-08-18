% Force v. Displacement Plot
close all; clc; clear all;

fileList = dir(fullfile(pwd,'*.mat'));
Fdimless = zeros(length(fileList),1);
gaps = zeros(length(fileList),1);
for i = 1:length(fileList)
    plotData_i = load(fileList(i).name);
    PlotData.Fdimless{i} = plotData_i.Fdim_save_mg;
    PlotData.gaps{i} = plotData_i.gap;
    PlotData.PE{i} = plotData_i.U_save_mg;
end

figure;
box on; hold on; grid on;


% scatter(PlotData.gaps{1}, PlotData.Fdimless{1}, 50, 'red', '^', 'filled');
scatter(1392, PlotData.Fdimless{2}, 50, 'blue', '^', 'filled');
scatter(5046, PlotData.Fdimless{4}, 50, 'r', '^', 'filled');
scatter(10216, PlotData.Fdimless{1}, 50, 'm', '^', 'filled');
scatter(28604, PlotData.Fdimless{3}, 50, 'c', '^', 'filled');




xlabel('{\boldmath $N$}','fontsize',20,'fontname','Times New Roman','interpreter','latex');
ylabel('{\boldmath $F/|F_0|$}','fontsize',24,'fontname','Times New Roman','interpreter','latex');
set(gcf,'position',[300,100,800,600]);
set(gca, 'LineWidth', 2.0 );
% set(gca, 'fontsize', 20.0 );
set(gca, 'XMinorTick', 'on');
% set(gca, 'Ticklength', [0.02;0.01] );
set(gca, 'YMinorTick', 'on');
% set(gca, 'xscale','log');
% set(gca, 'Ticklength', [0.02;0.01] );
% set(gca, 'XTick', 0:0.1:0.5 );
% set(gca, 'YTick', -0.4:.1:0.3 );

 ylim([0 1]);

legStr = {'$1392$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$5046$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$10216$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$28604$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
%     '$1392$ patches, $A = 0.05$, mountain-valley, $\sigma_1 = 1, \sigma_2 = -1$', ...
%     '$1562$ patches, $A = 0.1$, mountain-valley, $\sigma_1 = 1, \sigma_2 = -1$', ...
%     '$1704$ patches, $A = 0.15$, mountain-valley, $\sigma_1 = 1, \sigma_2 = -1$', ...
%     '$2732$ patches, $A = 0.20$, mountain-valley, $\sigma_1 = 1, \sigma_2 = -1$', ...
    };
% legStr = '1000 patches - Pos1Pos1';
legend(legStr,'location','northeast','fontsize',12,'interpreter','latex');



%{

% Raw Data - hardcoded. Should load from files instead
gap_2088_mv =  [10.0000
    9.0000
    8.0000
    7.0000
    6.0000
    5.0000
    4.0000
    3.0000
    2.0000
    1.0000
    0.9000
    0.8000
    0.7000
    0.6000
    0.5000
    0.4000
    0.3000
    0.2500];

gap_2088_tiptip = [10.0000
    9.0000
    8.0000
    7.0000
    6.0000
    5.0000
    4.0000
    3.0000
    2.0000
    1.0000
    0.9000
    0.8000
    0.7000
    0.6000
    0.5000
    0.4000
    0.3000
    0.2500];


gap_536_mv = [10.0000
    9.0000
    8.0000
    7.0000
    6.0000
    5.0000
    4.0000
    3.0000
    2.0000
    1.0000
    0.9000
    0.8000
    0.7000
    0.6000
    0.5000
    0.4000
    0.3000
    0.2000
    0.1000
    0.0900
    0.0800
    0.0700
    0.0600
    0.0500
    0.0400
    0.0300];

   gap_536_tiptip = [10.0000
    9.0000
    8.0000
    7.0000
    6.0000
    5.0000
    4.0000
    3.0000
    2.0000
    1.0000
    0.9000
    0.8000
    0.7000
    0.6000
    0.5000
    0.4000
    0.3000
    0.2000];

U_save_2088_mv_mg = [   -0.5366
   -0.5854
   -0.6440
   -0.7156
   -0.8051
   -0.9203
   -1.0742
   -1.2902
   -1.6166
   -2.1725
   -2.2510
   -2.3359
   -2.4280
   -2.5284
   -2.6385
   -2.7600
   -2.8950
   -2.9683];

U_save_2088_tiptip_mg = [   -0.5367
   -0.5855
   -0.6441
   -0.7157
   -0.8054
   -0.9207
   -1.0748
   -1.2913
   -1.6187
   -2.1775
   -2.2566
   -2.3421
   -2.4350
   -2.5364
   -2.6477
   -2.7709
   -2.9088
   -2.9797];

U_save_536_mv_mg = [   -0.3555
   -0.3878
   -0.4266
   -0.4741
   -0.5334
   -0.6097
   -0.7116
   -0.8548
   -1.0716
   -1.4441
   -1.4973
   -1.5550
   -1.6179
   -1.6869
   -1.7630
   -1.8478
   -1.9428
   -2.0506
   -2.1730
   -2.1860
   -2.1990
   -2.2121
   -2.2253
   -2.2386
   -2.2521
   -2.2660];

U_save_536_tiptip_mg = [-0.3557
   -0.3881
   -0.4270
   -0.4745
   -0.5340
   -0.6107
   -0.7132
   -0.8576
   -1.0771
   -1.4577
   -1.5125
   -1.5722
   -1.6375
   -1.7094
   -1.7894
   -1.8796
   -1.9832
   -2.1033];


Fdim_save_2088_mv_mg = [-0.0281
   -0.0334
   -0.0405
   -0.0500
   -0.0633
   -0.0828
   -0.1129
   -0.1632
   -0.2577
   -0.4739
   -0.5108
   -0.5525
   -0.6003
   -0.6554
   -0.7198
   -0.7956
   -0.8849
   -0.9344];


Fdim_save_2088_tiptip_mg = [   -0.0281
   -0.0335
   -0.0405
   -0.0500
   -0.0634
   -0.0829
   -0.1131
   -0.1636
   -0.2587
   -0.4772
   -0.5145
   -0.5569
   -0.6055
   -0.6616
   -0.7276
   -0.8070
   -0.9109
   -1.2164];

Fdim_save_536_mv_mg = [   -0.0281
   -0.0334
   -0.0405
   -0.0500
   -0.0633
   -0.0828
   -0.1130
   -0.1637
   -0.2598
   -0.4874
   -0.5275
   -0.5736
   -0.6269
   -0.6895
   -0.7637
   -0.8530
   -0.9622
   -1.0977
   -1.2763
   -1.2986
   -1.3221
   -1.3471
   -1.3733
   -1.4007
   -1.4289
   -1.4570];

Fdim_save_536_tiptip_mg = [   -0.0281
   -0.0335
   -0.0406
   -0.0501
   -0.0635
   -0.0832
   -0.1138
   -0.1653
   -0.2639
   -0.5021
   -0.5449
   -0.5945
   -0.6526
   -0.7220
   -0.8068
   -0.9141
   -1.0580
   -1.2888];


%}