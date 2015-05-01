close all;
clear all;
clc;

% Choose the number of basis functions
nBasis = 6;
averageError = 0.001;

%% Load data

% Wavelength sampling vector.
load('./Data/wave.mat');
% Camera matrix which represents the joint effects of the camera quantum
% efficiency and the external illumination.
load('./Data/camera.mat');
% Macbeth reflectance basis functions.
load('./Data/reflBasis.mat');
% Sample image data corrected for camera offset (dark voltage) and
% non-uniform illumination. Pixel values are therefore linear functions of
% the reflectance (which we are estimating).
load('./Data/macbethImage.mat');

reflBasis = reflBasis(:,1:nBasis);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

camera = camera*deltaL;

figure;
imshow(Img(:,:,5),[]);
title('Sample image channel');

h = size(Img,1);
w = size(Img,2);
nChannels = size(Img,3);


%% Solve spectral smoothness constraint only
%  This is a least-squares problem that can be solved quite efficiently
%  with Matlab functions.
alpha = 0.1;


tic
XlsCvx = spectralEstCVX(Img, camera, reflBasis, alpha, 0, 0);
tlsCvx = toc;


tic;
Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
measVals = reshape(Img,h*w,nChannels)';

mhat = [measVals; zeros(nWaves-1,size(measVals,2))];
Aalpha = [camera*reflBasis; sqrt(alpha)*Rlambda*reflBasis];
Xls = Aalpha\mhat;
Xls = reshape(Xls',[h w nBasis]);
tls = toc;



%% Spectral smoothness + box constraint
alpha = 0.1;
beta = 0;
gamma = 1;

tic
XBoxCvx = spectralEstCVX(Img, camera, reflBasis, alpha, beta, gamma);
tBoxCvx = toc;

tstart = tic;
[ XBoxAdmm, histBox ] = spectralEstADMM( Img, camera, reflBasis, alpha, beta, gamma,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',XBoxCvx,...
    'tolConv',averageError*(h*w));
tBoxAdmm = toc(tstart);

%% Spectral and spatial smoothness
alpha = 0.1;
beta = 0.1;
gamma = 0;

tic
XSpatialCvx = spectralEstCVX(Img, camera, reflBasis, alpha, beta, gamma);
tSpatialCvx = toc;


tstart = tic;
[ XSpatialAdmm, histSpatial ] = spectralEstADMM( Img, camera, reflBasis, alpha, beta, gamma,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',XSpatialCvx,...
    'tolConv',averageError*(h*w));
tSpatialAdmm = toc(tstart);

%% Spectral and spatial smoothness + box
alpha = 0.1;
beta = 0.1;
gamma = 1;

tic
XBoxSpatialCvx = spectralEstCVX(Img, camera, reflBasis, alpha, beta, gamma);
tBoxSpatialCvx = toc;

tstart = tic;
[ XBoxSpatialAdmm, histBoxSpatial ] = spectralEstADMM( Img, camera, reflBasis, alpha, beta, gamma,...
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

save('./Results/macbethResults.mat','tBoxAdmm','histBox','tSpatialAdmm','histSpatial','tBoxSpatialAdmm','histBoxSpatial',...
                'tBoxCvx','tSpatialCvx','tBoxSpatialCvx','tlsCvx','tls');
