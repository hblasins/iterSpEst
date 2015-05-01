close all;
clear variables;
clc;

%% This is a test script to show that the ADMM and cvx solutions are the same
%  For this purpose we create random camera, basis function and measurement
%  matrices.

% Choose problem size
nChannels = 7;
nWaves = 30;
nBasis = 5;
h = 5;
w = 5;

% Generate data
cameraMat = rand(nChannels,nWaves);
basisFcns = rand(nWaves,nBasis);
measImg = randn([h,w,nChannels]);


% Select solution regularizers. We only want to enforce spatial and spectral
% smoothness
alpha = 0.1;
beta = 0.1; 
gamma = 0;


%% CVX solution

Xcvx = spectralEstCVX(measImg, cameraMat, basisFcns, alpha, beta, gamma);

refl = basisFcns*Xcvx;

% Reflectance estimates
figure;
hold on; grid on; box on;
plot(refl);
xlabel('Wavelength, au');
title('Reflecctance estimates');

% Spatial arrangement of basis weights
figure;
for i=1:nBasis

    subplot(2,3,i);
    imagesc(reshape(Xcvx(i,:,:),h,w),[min(Xcvx(:)) max(Xcvx(:))]); axis image;
    title(sprintf('Basis %i',i));
    
end


%% ADMM solution

[Xadmm, hist]  = spectralEstADMM( measImg, cameraMat, basisFcns, alpha, beta, gamma,...
    'rescaleRho',true,...
    'maxIter',1000,...
    'verbose',true,...
    'reference',Xcvx);

% Convergence
figure; 
hold on; grid on; box on;
plot([hist.prRes hist.dualRes hist.conv]);
xlabel('Iteration');
legend('Primal','Dual','Solution');
title('ADMM residuals');

figure;
hold on; grid on; box on;
plot(hist.rho);
xlabel('Iteration');
title('Rho');

% Accuracy comparison with cvx
figure;
grid on; box on; hold on; axis square;
plot(Xadmm(:),Xcvx(:),'.');
xlabel('ADMM');
ylabel('CVX');
title(sprintf('ADMM vs. CVX accuracy (RMSE %f)',rms(Xcvx(:) - Xadmm(:))));



