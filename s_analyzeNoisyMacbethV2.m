% This is a demo script that runs the spectral estimatino algorithm on a
% set of images of a Macbeth test chart, acquied with different levels of
% read noise.
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

% Choose the number of basis functions
nBasis = 6;
wave = 380:5:950; wave = wave(:);

deltaL = wave(2) - wave(1);
nWaves = length(wave);


% The camera is the product of the quantum efficiency and the illuminant
fName = fullfile(iterSpEstRoot,'Parameters','sensorQe');
sensorQe = readSpectralQuantity(fName,wave);

fName = fullfile(iterSpEstRoot,'Parameters','illuminant');
illuminantEnergy = readSpectralQuantity(fName,wave);
illuminantQuanta = energy2Quanta(illuminantEnergy,wave);

% For numerical stability we normalize by scaleFac
camera = illuminantQuanta'*diag(sensorQe)*deltaL;


% Reference reflectance
fName = fullfile(iterSpEstRoot,'Parameters','macbethRefl');
reflRef = readSpectralQuantity(fName,wave);

% Reflectance basis functions
fName = fullfile(iterSpEstRoot,'Parameters','reflBasis');
reflBasis = readSpectralQuantity(fName,wave);
reflBasis = reflBasis(:,1:nBasis);



%% Make a prediction of the sensor responsivity on an ideally white,
%  uniform surface.
maxIter = 50;
alpha = 0.1;
beta = 0.1;
grType = 'anisotropic';

nPoints = 12;

RMSELs = zeros(nPoints,1);
RMSEBox = zeros(nPoints,1);
RMSESpatial = zeros(nPoints,1);
RMSEBoxSpatial = zeros(nPoints,1);

SNR = zeros(nPoints,1);

Images = cell(nPoints,1);
cameraMats = cell(nPoints,1);
for i=1:nPoints
    
    fName = fullfile(iterSpEstRoot,'Data',sprintf('noisyMacbethV2_%i',i));
    load(fName);
    
    h = size(Img,1);
    w = size(Img,2);
    Img = max(Img - repmat(shiftdim(cameraOffset,-2),[h w 1]),0);
    
    Images{i} = Img;
    cameraMats{i} = diag(cameraGain)*camera;
    
    SNR(i) = snr;

end
   
reflRefImg = reshape(reflRef',[4 6 nWaves]);
reflRefImg = imresize(reflRefImg,[h w],'nearest');
reflRefImg = subsampleImage(reflRefImg,4,6);

matlabpool open local
parfor i=1:nPoints

    Img = Images{i};
    cameraMat = cameraMats{i};
    

    % Least-squares
    [Xls] = spectralEstADMM(Img, cameraMat, reflBasis, alpha, 0, 0,...
        'maxIter',maxIter,...
        'verbose',true,...
        'gradient',grType);
    
    reflLs = reflBasis * reshape(Xls,h*w,nBasis)';
    reflLsImg = reshape(reflLs',[h,w,nWaves]);
    
    % For every pixel compute the l2 norm of the reflectance estimate error
    SEls = sum((subsampleImage(reflLsImg,4,6) - reflRefImg).^2,3);
    RMSELs(i) = sqrt(mean(SEls(:)));
    
    
    % showResults(subsampleImage(Xls,4,6), reflBasis, cameraMat, subsampleImage(Img,4,6), wave, reflRefImg);

    % Box
    [XBox] = spectralEstADMM(Img, cameraMat, reflBasis, alpha, 0, 1,...
        'maxIter',maxIter,...
        'verbose',true,...
        'gradient',grType);
    
    reflBox = reflBasis*reshape(XBox,h*w,nBasis)';
    reflBoxImg = reshape(reflBox',[h,w,nWaves]);
    
    SEBox = sum((subsampleImage(reflBoxImg,4,6) - reflRefImg).^2,3);
    RMSEBox(i) = sqrt(mean(SEBox(:)));
    
    
    % showResults(subsampleImage(XBox,4,6), reflBasis, cameraMat, subsampleImage(Img,4,6), wave, reflRefImg);

    % Spatial
    [XSpatial] = spectralEstADMM(Img, cameraMat, reflBasis, alpha, beta, 0,...
        'maxIter',maxIter,...
        'verbose',true,...
        'gradient',grType);
    
    reflSpatial = reflBasis*reshape(XSpatial,h*w,nBasis)';
    reflSpatialImg = reshape(reflSpatial',[h,w,nWaves]);
    
    SESpatial = sum((subsampleImage(reflSpatialImg,4,6) - reflRefImg).^2,3);
    RMSESpatial(i) = sqrt(mean(SESpatial(:)));
    
    
    % showResults(subsampleImage(XSpatial,4,6), reflBasis, cameraMat, subsampleImage(Img,4,6), wave, reflRefImg);

    % Box Spatial
    [XBoxSpatial] = spectralEstADMM(Img, cameraMat, reflBasis, alpha, beta, 1,...
        'maxIter',maxIter,...
        'verbose',true,...
        'gradient',grType);
    
    reflBoxSpatial = reflBasis*reshape(XBoxSpatial,h*w,nBasis)';
    reflBoxSpatialImg = reshape(reflBoxSpatial',[h,w,nWaves]);
    
    SEBoxSpatial = sum((subsampleImage(reflBoxSpatialImg,4,6) - reflRefImg).^2,3);
    RMSEBoxSpatial(i) = sqrt(mean(SEBoxSpatial(:)));
    
    % showResults(subsampleImage(XBoxSpatial,4,6), reflBasis, cameraMat, subsampleImage(Img,4,6), wave, reflRefImg);

    
end
matlabpool close

dName = fullfile(iterSpEstRoot,'Results');
if ~exist(dName,'dir');
    mkdir(dName);
end

fName = fullfile(iterSpEstRoot,'Results','noisyMacbethV2.mat');
save(fName,'SNR','RMSELs','RMSEBox','RMSESpatial','RMSEBoxSpatial');




