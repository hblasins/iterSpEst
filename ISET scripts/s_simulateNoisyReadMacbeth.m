% This is an ISET simulation script that generates a number of Macbeth
% chart scenes and simulates acquisition with cameras with different 
% levels of read noise.
%
% Copyright, Henryk Blasinski, 2015

close all;
clear variables;
clc;

s_initISET;
%%
wave = 380:5:950;

nWaves = length(wave);
deltaL = wave(2) - wave(1);


% Load the illuminant 
fName = fullfile(iterSpEstRoot,'Parameters','illuminant');
illuminantEnergy = ieReadSpectra(fName,wave);
illuminantQuanta = Energy2Quanta(wave,illuminantEnergy);
nChannels = size(illuminantEnergy,2);

% The camera is the product of the quantum efficiency and the illuminant
fName = fullfile(iterSpEstRoot,'Parameters','sensorQe');
sensorQe = readSpectralQuantity(fName,wave);

% Load the reference reflectances
fileName = fullfile(iterSpEstRoot,'Parameters','macbethRefl');
macbethRefl = ieReadSpectra(fileName,wave);


%% Simulate image acquisition with varying amount of read noise
cntr = 1;
readNoiseMult = logspace(0,2,20);

resFig = zeros(8,1);

for readNoise = [readNoiseMult(1) readNoiseMult(end-6:end)]

    Img = [];
    signal = [];
    signalWithNoise = [];
    
    cameraGain = zeros(nChannels,1);
    cameraOffset = zeros(nChannels,1);

    fprintf('Noise multiplier %.2f, channel: ', readNoise);
    
    for ch = 1:nChannels
    
        fprintf(' %i,',ch);
  
        [sensor,optics] = createFlea3Camera(wave);
    
    
        % Define a scene.
        scene = sceneCreate();
        scene = initDefaultSpectrum(scene,'custom',wave);
        macbethChartObject = macbethChartCreate(50,(1:24),sceneGet(scene,'spectrum'),'smallMacbethChart.mat');
        ill = illuminantCreate('equalphotons',wave);
        ill = illuminantSet(ill,'energy',illuminantEnergy(:,ch));
        scene = sceneCreateMacbeth(macbethChartObject,ill,scene);
     
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
        sensor = sensorSet(sensor,'noiseFlag',2);
        sensor = sensorSet(sensor,'Name',sprintf('Channel %i',ch));
    
        % Change sensor's noise characteristic
        pixel = sensorGet(sensor,'pixel');
        pixel = pixelSet(pixel,'readnoisevolts',pixelGet(pixel,'readnoisevolts')*readNoise);
        sensor = sensorSet(sensor,'pixel',pixel);
    
        % Compute and set the exposure duration
        cameraExposure = autoExposure(oi,sensor,0.9,'luminance');    
        sensor = sensorSet(sensor,'exposureTime',cameraExposure);
        
        [cameraGain(ch), cameraOffset(ch)] = sensorGetTotalGain(scene,oi,sensor);
    
        sensor = sensorCompute(sensor,oi);
        vcAddObject(sensor);
    
        image = double(sensorGet(sensor,'dv')/2^sensorGet(sensor,'nbits'));
        scaleFac = max(image(:));
        image = image/scaleFac;
        cameraGain(ch) = cameraGain(ch)/scaleFac;
        cameraOffset(ch) = cameraOffset(ch)/scaleFac;
        
        Img = cat(3,Img,image);
    
        % This will be necessary for input SNR calculations
        sensorNoiseless = sensorSet(sensor,'noiseFlag',0);
        sensorNoiseless = sensorCompute(sensorNoiseless,oi);
    
        signalWithNoise = cat(3,signalWithNoise,double(sensorGet(sensor,'volts')));
        signal = cat(3,signal,double(sensorGet(sensorNoiseless,'volts')));
        
        % Check if the measured data is linar with the prediction
        if cntr >= 1
            sz = sensorGet(sensor,'size');
            
            cp = [1 sz(1);
                sz(2) sz(1);
                sz(2) 1
                1 1];
            
            mVals = chartSelect(sensor,0,0,4,6,cp);
            mVals = mVals(:,1)/2^sensorGet(sensor,'nbits');
            mVals = mVals/scaleFac;
            
            pred = cameraGain(ch)*illuminantQuanta(:,ch)'*diag(sensorQe)*deltaL*macbethRefl + cameraOffset(ch);
            
            if resFig(cntr) == 0
                resFig(cntr) = figure;
                hold on; grid on; box on;
                xlabel('Simulated');
                ylabel('Modeled');
                title(sprintf('Noise multiplier %f',readNoise));
            else
                figure(resFig(cntr));
            end
            
            plot(mVals(:),pred(:),'.');
            drawnow;
        end
    
    end

    fprintf('\n');
    
    noise = signalWithNoise - signal;
    snr = 20*log10(mean(signal(:))/std(noise(:)));
    
    fName = fullfile(iterSpEstRoot,'Data',sprintf('noisyReadMacbeth_%i.mat',cntr));
    save(fName,'Img','cameraGain','cameraOffset','snr');
    cntr = cntr + 1;

end

