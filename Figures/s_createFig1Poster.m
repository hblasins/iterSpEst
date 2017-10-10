%% Run this script to generate Fig. 1 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015

wave = 400:10:700;

lw = 2;
fs = 10;


fName = fullfile(iterSpEstRoot,'Parameters','illuminant.mat');
light = ieReadSpectra(fName,wave);


figure; 
plot(wave,light(:,3),'lineWidth',lw);
% xlabel('Wavelength, nm','interpreter','latex','fontsize',fs);
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 5]);
set(gca,'YTick',0:0.01:0.03);
set(gca,'XTick',400:100:700);
plotTickLatex2D('yLabelDx',-0.02,'xLabelDy',0.01,'yDx',-1,'fontsize',fs-2);


fName = fullfile('~','Dropbox','MsVideo','Notes','ICIP15','Figures','LightsPoster.eps');
print('-depsc',fName);

fName = fullfile(iterSpEstRoot,'Parameters','sensorQe.mat');
qe = ieReadSpectra(fName,wave);

figure;
plot(wave,qe,'k','lineWidth',lw);
% xlabel('Wavelength, nm','interpreter','latex','fontsize',fs);
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 5]);
ylim([0 1.1*max(qe)]);
set(gca,'XTick',400:100:700);
plotTickLatex2D('yLabelDx',-0.02,'xLabelDy',0.01,'yDx',-1,'fontsize',fs-2);
fName = fullfile('~','Dropbox','MsVideo','Notes','ICIP15','Figures','CameraPoster.eps');
print('-depsc',fName);
