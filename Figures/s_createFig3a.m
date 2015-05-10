%% Run this script to generate Fig. 3 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

fName = fullfile(iterSpEstRoot,'Results','noisyScene.mat');
if ~exist(fName,'file');
    sName = fullfile(iterSpEstRoot,'s_analyzeNoisyScene.m');
    run(sName);
else
    load(fName);
end

fName = fullfile(iterSpEstRoot,'Data','sceneReference');
load(fName);

img = cat(3,sum(reflRefImg(:,:,41:end),3),sum(reflRefImg(:,:,21:40),3),sum(reflRefImg(:,:,1:20),3));
img = img/max(img(:));

%%

lw = 1.5;
fs = 6;
ms = 5;
f1 = figure;
hold on; grid on; box on;
plot(SNR,RMSELs,'g^-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSEBox,'rd-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSESpatial,'bs-','lineWidth',lw,'markersize',ms);
plot(SNR,RMSEBoxSpatial,'co-','lineWidth',lw,'markersize',ms);
xlabel('SNR, dB','fontsize',fs,'interpreter','latex');
ylabel('RMSE','fontsize',fs,'interpreter','latex');
ylim([0.4 1.4]);
set(gca','XTick',12:4:36);
set(gca','YTick',0.4:0.2:1.2);

sz = get(gcf,'Paperposition');
sz = 0.5*sz;
sz(4) = 0.5*sz(4);
set(gcf,'PaperPosition',sz);
plotTickLatex2D('yLabelDx',-0.02,'xLabelDy',0.01,'fontsize',fs-2);

% Add the image
ax = axes('Position',[0.5 0.4 0.5 0.5]);
imagesc(img);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis image;


fName = fullfile(iterSpEstRoot,'Figures','noisyScene.eps');
print('-depsc',fName);
