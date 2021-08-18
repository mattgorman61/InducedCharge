% Force v. Displacement Plot
close all; clc; clear all;

path1 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Tip-Tip\Pos1Neg1'; addpath(path1);
path2 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Tip-Tip\Pos1Pos1'; addpath(path2);
path3 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Tip-Tip\Pos1Pos10'; addpath(path3);
path4 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Mount-Val\Pos1Neg1'; addpath(path4);
path5 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Mount-Val\Pos1Pos1'; addpath(path5);
path6 = 'C:\Users\Matt Gorman\OneDrive\Documents\MATLAB\InducedCharge_git\3. Patch Resolution Study\Force Results\Mount-Val\Pos1Pos10'; addpath(path6);


fileList1 = dir(fullfile(path1,'*.mat'));
for i = 1:length(fileList1)
    plotData1_i = load(fileList1(i).name);
    PlotData1.Fdimless{i} = plotData1_i.Fdim_save_mg;
    PlotData1.gaps{i} = plotData1_i.gap;
    PlotData1.PE{i} = plotData1_i.U_save_mg;
end

fileList2 = dir(fullfile(path2,'*.mat'));
for i = 1:length(fileList2)
    plotData2_i = load(fileList2(i).name);
    PlotData2.Fdimless{i} = plotData2_i.Fdim_save_mg;
    PlotData2.gaps{i} = plotData2_i.gap;
    PlotData2.PE{i} = plotData2_i.U_save_mg;
end

fileList3 = dir(fullfile(path3,'*.mat'));
for i = 1:length(fileList3)
    plotData3_i = load(fileList3(i).name);
    PlotData3.Fdimless{i} = plotData3_i.Fdim_save_mg;
    PlotData3.gaps{i} = plotData3_i.gap;
    PlotData3.PE{i} = plotData3_i.U_save_mg;
end

fileList4 = dir(fullfile(path4,'*.mat'));
for i = 1:length(fileList4)
    plotData4_i = load(fileList4(i).name);
    PlotData4.Fdimless{i} = plotData4_i.Fdim_save_mg;
    PlotData4.gaps{i} = plotData4_i.gap;
    PlotData4.PE{i} = plotData4_i.U_save_mg;
end

fileList5 = dir(fullfile(path5,'*.mat'));
for i = 1:length(fileList5)
    plotData5_i = load(fileList5(i).name);
    PlotData5.Fdimless{i} = plotData5_i.Fdim_save_mg;
    PlotData5.gaps{i} = plotData5_i.gap;
    PlotData5.PE{i} = plotData5_i.U_save_mg;
end

fileList6 = dir(fullfile(path6,'*.mat'));
for i = 1:length(fileList6)
    plotData6_i = load(fileList6(i).name);
    PlotData6.Fdimless{i} = plotData6_i.Fdim_save_mg;
    PlotData6.gaps{i} = plotData6_i.gap;
    PlotData6.PE{i} = plotData6_i.U_save_mg;
end

figure;
box on; hold on; grid on;


scatter(1392, PlotData1.Fdimless{2}, 100, 'r', '^', 'filled');
scatter(5046, PlotData1.Fdimless{4}, 100, 'r', '^', 'filled');
scatter(10216, PlotData1.Fdimless{1}, 100, 'r', '^', 'filled');
scatter(28604, PlotData1.Fdimless{3}, 100, 'r', '^', 'filled');

scatter(1392, PlotData2.Fdimless{2}, 100, 'b', '^', 'filled');
scatter(5046, PlotData2.Fdimless{4}, 100, 'b', '^', 'filled');
scatter(10216, PlotData2.Fdimless{1}, 100, 'b', '^', 'filled');
scatter(28604, PlotData2.Fdimless{3}, 100, 'b', '^', 'filled');

scatter(1392, PlotData3.Fdimless{2}, 100, 'k', '^', 'filled');
scatter(5046, PlotData3.Fdimless{4}, 100, 'k', '^', 'filled');
scatter(10216, PlotData3.Fdimless{1}, 100, 'k', '^', 'filled');
scatter(28604, PlotData3.Fdimless{3}, 100, 'k', '^', 'filled');

scatter(1392, PlotData4.Fdimless{2}, 100, 'r', 'x');
scatter(5046, PlotData4.Fdimless{4}, 100, 'r', 'x');
scatter(10216, PlotData4.Fdimless{1}, 100, 'r', 'x');
scatter(28604, PlotData4.Fdimless{3}, 100, 'r', 'x');

scatter(1392, PlotData5.Fdimless{2}, 100, 'b', 'x');
scatter(5046, PlotData5.Fdimless{4}, 100, 'b', 'x');
scatter(10216, PlotData5.Fdimless{1}, 100, 'b', 'x');
scatter(28604, PlotData5.Fdimless{3}, 100, 'b', 'x');

scatter(1392, PlotData6.Fdimless{2}, 100, 'k', 'x');
scatter(5046, PlotData6.Fdimless{4}, 100, 'k', 'x');
scatter(10216, PlotData6.Fdimless{1}, 100, 'k', 'x');
scatter(28604, PlotData6.Fdimless{3}, 100, 'k', 'x');

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

%  ylim([-2 0]);
xlim([0,30000]);

legStr1 = {'$1392$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$5046$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$10216$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$28604$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = -1$', ...
};

legStr2 = {'$1392$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$5046$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$10216$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$28604$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 1$', ...
};

legStr3 = {'$1392$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$5046$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$10216$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$28604$ patches, $A = 0.05$, tip-tip, $\sigma_1 = 1, \sigma_2 = 10$', ...
};

legStr4 = {'$1392$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$5046$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$10216$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = -1$', ...
    '$28604$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = -1$', ...
};

legStr5 = {'$1392$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$5046$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$10216$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
    '$28604$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 1$', ...
};

legStr6 = {'$1392$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$5046$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$10216$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 10$', ...
    '$28604$ patches, $A = 0.05$, m-v, $\sigma_1 = 1, \sigma_2 = 10$', ...
};

legend([legStr1,legStr2,legStr3,legStr4,legStr5,legStr6],'location','northeast','fontsize',12,'interpreter','latex','numColumns',2);

