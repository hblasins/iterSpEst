%% Run this script to generate Fig. 3 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

fName = fullfile(iterSpEstRoot,'Results','noisyMacbeth.mat');
if ~exist(fName,'file');
    sName = fullfile(iterSpEstRoot,'s_analyzeNoisyMacbeth.m');
    run(sName);
else
    load(fName);
end


%%

lw = 1.5;
fs = 6;
ms = 5;
figure;
hold on; grid on; box on;
plot(SNR,RMSELs,'g^-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSEBox,'rd-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSESpatial,'bs-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSEBoxSpatial,'co-','lineWidth',lw,'markersize',ms);
xlabel('SNR, dB','fontsize',fs,'interpreter','latex');
ylabel('RMSE','fontsize',fs,'interpreter','latex');
ylim([0.2 1.4]);
set(gca','XTick',12:4:36);
set(gca','YTick',0.4:0.2:1.2);

lg = legend({'Least-squares','Box','Spatial','Box Spatial'},...
    'location','northeast','interpreter','latex','fontsize',fs-2);
lgSz = get(lg,'Position');
lgSz(3) = lgSz(3) *1.3;
lgSz(4) = 0.6*lgSz(4);
lgSz(1) = lgSz(1) - 0.15;
lgSz(2) = lgSz(2) - 0.15;
set(lg,'Position',lgSz);

sz = get(gcf,'Paperposition');
sz = 0.5*sz;
sz(4) = 0.5*sz(4);
set(gcf,'PaperPosition',sz);
plotTickLatex2D('yLabelDx',-0.02,'xLabelDy',0.01,'fontsize',fs-2);
fName = fullfile(iterSpEstRoot,'Figures','noisyMacbeth.eps');
print('-depsc',fName);
