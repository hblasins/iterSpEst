% This is an ISET simulation script that generates a number of natural
% scenes and simulates acquisition with cameras with different levels of 
% read noise.
%
% Copyright, Henryk Blasinski, 2015

close all;
clear variables;
clc;

s_initISET;
%%

nChannels = 7;
nBasis = 6;
wave = 400:5:700;

nWaves = length(wave);
deltaL = wave(2) - wave(1);


% Load the illuminant 
fName = fullfile(iterSpEstRoot,'Parameters','illuminant');
illuminantEnergy = ieReadSpectra(fName,wave);
illuminantEnergy = illuminantEnergy(:,1:nChannels);
illuminantQuanta = Energy2Quanta(wave,illuminantEnergy);


% The camera is the product of the quantum efficiency and the illuminant
fName = fullfile(iterSpEstRoot,'Parameters','sensorQe');
sensorQe = readSpectralQuantity(fName,wave);

% Load the reference reflectances
fileName = fullfile(iterSpEstRoot,'Parameters','macbethRefl');
macbethRefl = ieReadSpectra(fileName,wave);

% Load the basis functions
fileName = fullfile(iterSpEstRoot,'Parameters','reflBasis');
reflBasis = ieReadSpectra(fileName,wave);
reflBasis = reflBasis(:,1:nBasis);


%% Simulations
cntr = 1;
nPoints = 8;
cameraExposure = zeros(nPoints,1);

% Create a scene.
fName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
refScene = sceneFromFile(fName,'multispectral',[],[],wave);

for lightMult = logspace(0,-2,nPoints);

    Img = [];
    signal = [];
    signalWithNoise = [];

    cameraGain = zeros(nChannels,1);
    cameraOffset = zeros(nChannels,1);

    fprintf('Light multiplier %.3f, channel: ', lightMult);
    for ch = 1:nChannels
        
        fprintf(' %i,',ch);
        [sensor, optics] = createFlea3Camera(wave);
    
            
        
        scene = sceneAdjustIlluminant(refScene,lightMult*illuminantEnergy(:,ch),0); 
        scene = sceneSet(scene,'fov',1);
        scene = sceneSet(scene,'distance',1);
        scene = sceneSet(scene,'name',sprintf('Channel %i',ch));
        vcAddObject(scene);
    
        % Compute the optical image
        oi = oiCreate();
        oi = oiSet(oi,'optics',optics);
        oi = oiSet(oi,'Name',sprintf('Channel %i',ch));
        oi = oiCompute(scene,oi);
        vcAddObject(oi);
    
        % Compute the sensor image
        FOV = [sceneGet(scene,'fov horizontal') sceneGet(scene,'fov vertical')];
        sensor = sensorSetSizeToFOV(sensor,FOV,scene,oi);
        sensor = sensorSet(sensor,'Name',sprintf('Channel %i',ch));

        if cntr == 1
            cameraExposure(ch) = autoExposure(oi,sensor,0.9,'luminance');
        end
        sensor = sensorSet(sensor,'exposureTime',cameraExposure(ch));
        [cameraGain(ch), cameraOffset(ch)] = sensorGetTotalGain(scene,oi,sensor);

        sensor = sensorCompute(sensor,oi);
        vcAddObject(sensor);

        image = double(sensorGet(sensor,'dv'))/(2^sensorGet(sensor,'nbits'));
        scaleFac = max(image(:));
        image = image/scaleFac;
        cameraGain(ch) = cameraGain(ch)/scaleFac;
        cameraOffset(ch) = cameraOffset(ch)/scaleFac;
        
        ROI = image(23:36,218:236);
        
        Img = cat(3,Img,image);

        % Noise computation
        image = image - cameraOffset(ch);
        scaleFac = max(image(:));
        image = image/scaleFac;
        
        sensorNoiseless = sensorSet(sensor,'noise flag',0);
        sensorNoiseless = sensorCompute(sensorNoiseless,oi);

        [gainNoiseless, offsetNoiseless] = sensorGetTotalGain(scene,oi,sensorNoiseless);
        imageNoisless = sensorGet(sensorNoiseless,'dv')/2^sensorGet(sensor,'nbits') - offsetNoiseless;
        imageNoisless = imageNoisless/scaleFac;
        
        signal = cat(3,signal,imageNoisless);
        signalWithNoise = cat(3,signalWithNoise,image);
        
    
    end
    
    fprintf('\n');

    noise = signalWithNoise - signal;
    snr = 20*log10(mean(signal(:))/std(noise(:)));

    
    cameraGain = cameraGain*lightMult;
    
    fName = fullfile(iterSpEstRoot,'Data',sprintf('noisyLightScene_%i.mat',cntr)); 
    save(fName,'Img','cameraGain','cameraOffset','snr');
    cntr = cntr + 1;

end

%% Recover the full reflectance spectrum
% We need to simulate the entire camera pipeline to account for blurring 

reflRefImg = [];

fName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs.mat');
scene = sceneFromFile(fName,'multispectral',[],[],wave);
scene = sceneSet(scene,'fov',1);
scene = sceneSet(scene,'distance',1);
vcAddObject(scene);

ill = sceneGet(scene,'illuminant photons');

fprintf('Reference channel: ');
for ch = 1:nWaves
    
    fprintf(' %i,',ch);

    
    [sensor,optics] = createFlea3Camera(wave);
    
    pixel = sensorGet(sensor,'pixel');
    pixel = pixelSet(pixel,'readnoisevolts',pixelGet(pixel,'readnoisevolts'));
    qe = zeros(nWaves,1);
    qe(ch) = 1;
    pixel = pixelSet(pixel,'pixelspectralqe',qe);
    sensor = sensorSet(sensor,'pixel',pixel);
    
    % Compute the optical image
    oi = oiCreate();
    oi = oiSet(oi,'optics',optics);
    oi = oiSet(oi,'Name',sprintf('Channel %i',ch));
    oi = oiCompute(scene,oi);
    vcAddObject(oi);
    
    % Compute the sensor image
    FOV = [sceneGet(scene,'fov horizontal') sceneGet(scene,'fov vertical')];
    sensor = sensorSetSizeToFOV(sensor,FOV,scene,oi);
    sensor = sensorSet(sensor,'Name',sprintf('Channel %i',ch));
    sensor = sensorSet(sensor,'noiseFlag',-1);
    
    cameraExposure = autoExposure(oi,sensor,0.9,'luminance');    
    sensor = sensorSet(sensor,'exposureTime',cameraExposure);
    cameraGain = sensorGetTotalGain(scene,oi,sensor);
    
    sensor = sensorCompute(sensor,oi);
    vcAddObject(sensor);
    
    image = sensorGet(sensor,'volts')/pixelGet(pixel,'voltageSwing');
    
    % Divide by stuff to get reflectance
    image = image/deltaL/cameraGain/ill(ch);
    
    reflRefImg = cat(3,reflRefImg,image);
    
end

fprintf('\n');

fName = fullfile(iterSpEstRoot,'Data','sceneReference.mat'); 
save(fName,'reflRefImg');
