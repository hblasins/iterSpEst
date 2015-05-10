%% Run this script to generate Fig. 1 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

fName = fullfile(iterSpEstRoot,'Results','sampleResults.mat');
if ~exist(fName,'file');
    sName = fullfile(iterSpEstRoot,'s_analyzeAccuracy');
    run(sName);
else
    load(fName);
end


% Crop images
Xls = Xls(200:end,1:end-240,:);
XBoxAdmm = XBoxAdmm(200:end,1:end-240,:);
XSpatialAdmm = XSpatialAdmm(200:end,1:end-240,:);
XBoxSpatialAdmm = XBoxSpatialAdmm(200:end,1:end-240,:);

%% Spatial plots

x = [ 11 228 7 241 344 574];
y = [ 291 517 638 819 294 514] - 200;
scale = 1;

for wg=1:2
    maxVal = max(max([Xls(:,:,wg), XBoxAdmm(:,:,wg), XSpatialAdmm(:,:,wg) XBoxSpatialAdmm(:,:,wg)]));
    minVal = min(min([Xls(:,:,wg), XBoxAdmm(:,:,wg), XSpatialAdmm(:,:,wg) XBoxSpatialAdmm(:,:,wg)]));
    
    i1 = (Xls(:,:,wg) - minVal)/(maxVal-minVal);
    i2 = (XBoxAdmm(:,:,wg) - minVal)/(maxVal-minVal);
    i3 = (XSpatialAdmm(:,:,wg) - minVal)/(maxVal-minVal);
    i4 = (XBoxSpatialAdmm(:,:,wg) - minVal)/(maxVal-minVal);
    
    xp = x*scale;
    yp = y*scale;
    
    f1 = figure;
    hndl = imshow(imresize(i1,scale,'nearest'),'Border','tight','InitialMagnification',100);
    colormap(gray(256));
    hold on;
    plt = plot(300 - 400*i1(300,:),'Color',[0.95 0.95 0.95],'LineWidth',2);
    if wg == 1
        r = zeros(3,1);
        t = zeros(3,1);
        for j=1:2
            ind = 2*(j-1)+1;
            if j==2
                r(j) = rectangle('Position',[xp(ind) yp(ind) xp(ind+1)-xp(ind) yp(ind+1)-yp(ind)],'EdgeColor','black','LineWidth',4);
                t(j) = text(xp(ind)+30,yp(ind)+30,sprintf('#%i',j),'FontSize',28);
            else
                r(j) = rectangle('Position',[xp(ind) yp(ind) xp(ind+1)-xp(ind) yp(ind+1)-yp(ind)],'EdgeColor',[0.95 0.95 0.95],'LineWidth',4);
                t(j) = text(xp(ind)+30,yp(ind)+30,sprintf('#%i',j),'FontSize',28,'Color',[0.95 0.95 0.95]);
            end
        end
    end
    
    fName = fullfile(iterSpEstRoot,'Figures',sprintf('LSimg_w%i.eps',wg));
    print('-depsc',fName);
    
    colormap(gray(256));
    left=100; bottom=100 ; width=100 ; height=500;
    pos=[left bottom width height];
    axis off;
    set(hndl,'Visible','off');
    set(plt,'Visible','off');
    for j=1:3, set(r(j),'Visible','off'); set(t(j),'Visible','off'); end;
    labels = linspace(ceil(minVal),floor(maxVal),floor(maxVal) - ceil(minVal) + 1);
    dist = (labels - minVal)/(maxVal - minVal);
    sz = get(f1,'PaperPosition');
    sz(3) = 0.1*sz(3);
    set(f1,'PaperPosition',sz);
    set(f1,'OuterPosition',pos);
    cb = colorbar([0.1 0.1  0.5  0.8],'YTick',dist,'YTickLabel',num2str(labels'),'FontSize',16);

    fName = fullfile(iterSpEstRoot,'Figures',sprintf('Colorbar_w%i.eps',wg));
    print('-depsc',fName);
    
    
    figure;
    imshow(imresize(i2,scale,'nearest'),'Border','tight','InitialMagnification',100);
    colormap(gray(256));
    hold on;
    plot(300 - 400*i2(300,:),'Color',[0.95 0.95 0.95],'LineWidth',2);
    
    fName = fullfile(iterSpEstRoot,'Figures',sprintf('Boximg_w%i.eps',wg));
    print('-depsc',fName);
    
    figure;
    imshow(imresize(i3,scale,'nearest'),'Border','tight','InitialMagnification',100);
    colormap(gray(256));
    hold on;
    plot(300 - 400*i3(300,:),'Color',[0.95 0.95 0.95],'LineWidth',2);
    fName = fullfile(iterSpEstRoot,'Figures',sprintf('Spatialimg_w%i.eps',wg));
    print('-depsc',fName);
    
    figure;
    imshow(imresize(i4,scale,'nearest'),'Border','tight','InitialMagnification',100);
    colormap(gray(256));
    hold on;
    plot(300 - 400*i4(300,:),'Color',[0.95 0.95 0.95],'LineWidth',2);
    fName = fullfile(iterSpEstRoot,'Figures',sprintf('BoxSpatialimg_w%i.eps',wg));
    print('-depsc',fName);
end
