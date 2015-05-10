% This is a demo script that runs the spectral estimatino algorithm on a
% subsamples image of a Macbeth chart and compares the execution time of a
% conventional approach (CVX) and the proposed iterative method (ADMM).
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

% Choose the number of basis functions
nBasis = 6;
averageError = 0.001;
wave = 380:5:950; wave = wave(:);

deltaL = wave(2) - wave(1);
nWaves = length(wave);


% The camera is the product of the quantum efficiency and the illuminant
fName = fullfile(iterSpEstRoot,'Parameters','sensorQe');
sensorQe = readSpectralQuantity(fName,wave);

fName = fullfile(iterSpEstRoot,'Parameters','illuminant');
illuminantEnergy = readSpectralQuantity(fName,wave);
illuminantQuanta = energy2Quanta(illuminantEnergy,wave);

camera = illuminantQuanta'*diag(sensorQe)*deltaL;

% Reference reflectance
fName = fullfile(iterSpEstRoot,'Parameters','macbethRefl');
reflRef = readSpectralQuantity(fName,wave);

% Reflectance basis functions
fName = fullfile(iterSpEstRoot,'Parameters','reflBasis');
reflBasis = readSpectralQuantity(fName,wave);
reflBasis = reflBasis(:,1:nBasis);

% Load the sample macbeth Image
fName = fullfile(iterSpEstRoot,'Data','macbethImage');
load(fName);
h = size(Img,1);
w = size(Img,2);
nChannels = size(Img,3);

Img = max(Img - repmat(shiftdim(cameraOffset',-2),[h w 1]),0);
cameraMat = diag(cameraGain)*camera;

reflRefImg = imresize(reshape(reflRef',[4 6 nWaves]),[h w],'nearest');


figure;
imshow(Img(:,:,5),[]);
title('Sample image channel');




%% Solve spectral smoothness constraint only
%  This is a least-squares problem that can be solved quite efficiently
%  with Matlab functions.
alpha = 0.1;


tic
XlsCvx = spectralEstCVX(Img, cameraMat, reflBasis, alpha, 0, 0);
tlsCvx = toc;


tic;
Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
measVals = reshape(Img,h*w,nChannels)';

mhat = [measVals; zeros(nWaves-1,size(measVals,2))];
Aalpha = [cameraMat*reflBasis; sqrt(alpha)*Rlambda*reflBasis];
Xls = Aalpha\mhat;
Xls = reshape(Xls',[h w nBasis]);
tls = toc;



%% Spectral smoothness + box constraint
alpha = 0.1;
beta = 0;
gamma = 1;

tic
XBoxCvx = spectralEstCVX(Img, cameraMat, reflBasis, alpha, beta, gamma);
tBoxCvx = toc;

tstart = tic;
[ XBoxAdmm, histBox ] = spectralEstADMM( Img, cameraMat, reflBasis, alpha, beta, gamma,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',XBoxCvx,...
    'tolConv',averageError*(h*w));
tBoxAdmm = toc(tstart);

showResults( XBoxCvx, reflBasis, cameraMat, Img, wave, reflRefImg );

%% Spectral and spatial smoothness
alpha = 0.1;
beta = 0.1;
gamma = 0;

tic
XSpatialCvx = spectralEstCVX(Img, cameraMat, reflBasis, alpha, beta, gamma);
tSpatialCvx = toc;


tstart = tic;
[ XSpatialAdmm, histSpatial ] = spectralEstADMM( Img, cameraMat, reflBasis, alpha, beta, gamma,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',XSpatialCvx,...
    'tolConv',averageError*(h*w));
tSpatialAdmm = toc(tstart);

showResults( XSpatialAdmm, reflBasis, cameraMat, Img, wave, reflRefImg );


%% Spectral and spatial smoothness + box
alpha = 0.1;
beta = 0.1;
gamma = 1;

tic
XBoxSpatialCvx = spectralEstCVX(Img, cameraMat, reflBasis, alpha, beta, gamma);
tBoxSpatialCvx = toc;

tstart = tic;
[ XBoxSpatialAdmm, histBoxSpatial ] = spectralEstADMM( Img, cameraMat, reflBasis, alpha, beta, gamma,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',XSpatialCvx,...
    'tolConv',averageError*(h*w));
tBoxSpatialAdmm = toc(tstart);

% Timing results
fprintf('Spectrally smooth (least-squares):            %f sec.\n',tls);
fprintf('Spectrally smooth + box (ADMM):               %f sec\n',tBoxAdmm);
fprintf('Spectrally and spatially smooth (ADMM):       %f sec\n',tSpatialAdmm);
fprintf('Spectrally and spatially smooth + Box (ADMM): %f sec\n',tBoxSpatialAdmm);

dName = fullfile(iterSpEstRoot,'Results');
if ~exist(dName,'dir');
    mkdir(dName);
end

fName = fullfile(iterSpEstRoot,'Results','macbethResults.mat');
save(fName,'tBoxAdmm','histBox','tSpatialAdmm','histSpatial','tBoxSpatialAdmm','histBoxSpatial',...
                'tBoxCvx','tSpatialCvx','tBoxSpatialCvx','tlsCvx','tls');
