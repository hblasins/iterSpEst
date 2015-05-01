close all;
clear all;
clc;

if ~exist('./Results/sampleResults.mat','file');
    run('../s_analyzeAccuracy.m');
else
    load('./Results/sampleResults.mat');
end

if ~exist('./Figures','dir');
    mkdir('./Figures');
end

% Macbeth patch ID's that are analyzed
patchIDs = [3, 4, 7];
load('./Data/reflBasis.mat');
load('./Data/wave.mat');
load('./Data/macbethRefl.mat');
nWaves = length(wave);

% Crop images
Xls = Xls(200:end,1:end-240,:);
XBoxAdmm = XBoxAdmm(200:end,1:end-240,:);
XSpatialAdmm = XSpatialAdmm(200:end,1:end-240,:);
XBoxSpatialAdmm = XBoxSpatialAdmm(200:end,1:end-240,:);

[h, w, nBasis] = size(Xls);
reflBasis = reflBasis(:,1:nBasis);

XlsRefl = reshape((reflBasis*reshape(Xls,h*w,nBasis)')',[h w nWaves]);
XBoxAdmmRefl = reshape((reflBasis*reshape(XBoxAdmm,h*w,nBasis)')',[h w nWaves]);
XSpatialAdmmRefl = reshape((reflBasis*reshape(XSpatialAdmm,h*w,nBasis)')',[h w nWaves]);
XBoxSpatialAdmmRefl = reshape((reflBasis*reshape(XBoxSpatialAdmm,h*w,nBasis)')',[h w nWaves]);


%% Spectral plots

x = [ 11 228 7 241 344 574];
y = [ 291 517 638 819 294 514] - 200;
scale = 1;

%% Spectral curves
clear sz;
lw = 3;

for j=1:2
    indx = 2*(j-1)+1;
    subRefl = XlsRefl(y(indx):5:y(indx+1),x(indx):5:x(indx+1),:);
    subRefl = reshape(subRefl,size(subRefl,1)*size(subRefl,2),nWaves)';
    
    figure;
    hold on; grid on; box on;
    plot(wave,subRefl,'k');
    plot(wave,macbethRefl(:,patchIDs(j)),'r--','LineWidth',lw);
    plot(wave,zeros(size(wave)),'g-','LineWidth',lw);
    plot(wave,ones(size(wave)),'g-','LineWidth',lw);
    xlim([min(wave) max(wave)]);
    set(gca,'XTick',400:100:1000);
    set(gca,'YTick',0:0.5:1);

    ylim([-0.1 1.1]);
    
    if ~exist('sz','var');
        sz = get(gcf,'PaperPosition');
        sz(3:4) = sz(3:4)*0.5;
        sz(4) = 0.5*sz(4);
    end
    set(gcf,'PaperPosition',sz);
    plotTickLatex2D('xDy',-0.03);
    print('-depsc',sprintf('./Figures/RealLS_%i.eps',patchIDs(j)));
end


for j=1:2
    indx = 2*(j-1)+1;
    subRefl = XBoxAdmmRefl(y(indx):5:y(indx+1),x(indx):5:x(indx+1),:);
    subRefl = reshape(subRefl,size(subRefl,1)*size(subRefl,2),nWaves)';
    
    figure;
    hold on; grid on; box on;
    plot(wave,subRefl,'k');
    plot(wave,macbethRefl(:,patchIDs(j)),'r--','LineWidth',lw);
    plot(wave,zeros(size(wave)),'g-','LineWidth',lw);
    plot(wave,ones(size(wave)),'g-','LineWidth',lw);
    xlim([min(wave) max(wave)]);
    set(gca,'XTick',400:100:1000);
    set(gca,'YTick',0:0.5:1);

    ylim([-0.1 1.1]);
    
    set(gcf,'PaperPosition',sz);
    plotTickLatex2D('xDy',-0.03);
    print('-depsc',sprintf('./Figures/RealBox_%i.eps',patchIDs(j)));
end

for j=1:2
    indx = 2*(j-1)+1;
    subRefl = XSpatialAdmmRefl(y(indx):5:y(indx+1),x(indx):5:x(indx+1),:);
    subRefl = reshape(subRefl,size(subRefl,1)*size(subRefl,2),nWaves)';
    
    figure;
    hold on; grid on; box on;
    plot(wave,subRefl,'k');
    plot(wave,macbethRefl(:,patchIDs(j)),'r--','LineWidth',lw);
    plot(wave,zeros(size(wave)),'g-','LineWidth',lw);
    plot(wave,ones(size(wave)),'g-','LineWidth',lw);
    xlim([min(wave) max(wave)]);
    set(gca,'XTick',400:100:1000);   
    set(gca,'YTick',0:0.5:1);

    ylim([-0.1 1.1]);
    
    set(gcf,'PaperPosition',sz);
    plotTickLatex2D('xDy',-0.03);
    print('-depsc',sprintf('./Figures/RealSpatial_%i.eps',patchIDs(j)));
end

for j=1:2
    indx = 2*(j-1)+1;
    subRefl = XBoxSpatialAdmmRefl(y(indx):5:y(indx+1),x(indx):5:x(indx+1),:);
    subRefl = reshape(subRefl,size(subRefl,1)*size(subRefl,2),nWaves)';
    
    figure;
    hold on; grid on; box on;
    plot(wave,subRefl,'k');
    plot(wave,macbethRefl(:,patchIDs(j)),'r--','LineWidth',lw);
    plot(wave,zeros(size(wave)),'g-','LineWidth',lw);
    plot(wave,ones(size(wave)),'g-','LineWidth',lw);
    xlim([min(wave) max(wave)]);
    set(gca,'XTick',400:100:1000);
    set(gca,'YTick',0:0.5:1);

    ylim([-0.1 1.1]);
    
    set(gcf,'PaperPosition',sz);
    plotTickLatex2D('xDy',-0.03);
    print('-depsc',sprintf('./Figures/RealBoxSpatial_%i.eps',patchIDs(j)));
end

