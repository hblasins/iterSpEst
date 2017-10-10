%% Run this script to generate Fig. 3 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

fName = fullfile(iterSpEstRoot,'Results','noisyLightScene.mat');
if ~exist(fName,'file');
    sName = fullfile(iterSpEstRoot,'s_analyzeNoisyLightScene.m');
    run(sName);
else
    load(fName);
end

fName = fullfile(iterSpEstRoot,'Data','sceneReference');
load(fName);

img = cat(3,sum(reflRefImg(:,:,41:end),3),sum(reflRefImg(:,:,21:40),3),sum(reflRefImg(:,:,1:20),3));
img = (img/max(img(:))).^(1/2.2);

%%

lw = 5;
fs = 20;
ms = 10;
f1 = figure;
hold on; grid on; box on;
plot(SNR,RMSELs,'g^-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSEBox,'rd-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSESpatial,'bs-','lineWidth',lw,'markersize',2*ms);
plot(SNR,RMSEBoxSpatial,'co-','lineWidth',lw,'markersize',ms);
xlabel('Measurement SNR, dB','fontsize',fs,'interpreter','latex');
ylabel('Reflectance RMSE','fontsize',fs,'interpreter','latex');
ylim([0.4 1.4]);
set(gca,'XTick',12:2:36);
set(gca,'YTick',0.4:0.2:1.2);
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 24 15]);
plotTickLatex2D('yLabelDx',-0.02,'xLabelDy',0.01,'fontsize',fs-2);

% Add the image
ax = axes('Position',[0.5 0.5 0.4 0.4]);
imagesc(img);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
axis image;


fName = fullfile('~','Dropbox','MsVideo','Notes','ICIP15','Figures','noisyLightScenePoster.eps');
print('-depsc',fName);
