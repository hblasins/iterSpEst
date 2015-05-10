function [sensor, optics] = createFlea3Camera(wave)

% function [sensor, optics] = createFlea3Camera(wave)
%
% This function creates an ISET model of a PtGrey Flea3 camera. Most data
% is publicly available on various datasheets. The quantum efficiency was 
% measured.
%
% Copyright, Henryk Blasinski, 2015

nWaves = length(wave);

% Here are some of the key pixel properties of the PointGrey Flea3 camera
% with an ONSemi VITA 1300 1.3 Megapixel sensor.
wellCapacity   = 13700;                     % Electrons
conversiongain = 90*1e-6;                   % Volts/electron 
voltageSwing   = conversiongain*wellCapacity;
fillfactor     = 0.5;                       % This is a made up number
pixelSize      = 4.8*1e-6;                  % Meters
darkvoltage    = conversiongain*4.5;        % Volts/sec
readnoise      = 0.00096;                   % Volts
rows = 1024;                                 
cols = 1280; 
dsnu =  conversiongain*30;                  % Volts (dark signal non-uniformity)
prnu = 2;                                   % Percent (ranging between 0 and 100) photodetector response non-uniformity
quantizationMethod = '8bit';
analogGain   = 1/1.3089;                    % Gain and offset determined via calibration
analogOffset = 0.0255*1.3089*voltageSwing;               

fName = fullfile(iterSpEstRoot,'Parameters','sensorQe');
qe = readSpectralQuantity(fName,wave); % Spectral qe determined via calibration

% Some parameters of the lens
focalLength = 0.07;
fNumber = 4;



sensor = sensorCreate('Monochrome');
pixel =  sensorGet(sensor,'pixel');

% We set these properties here
pixel = pixelSet(pixel,'sizesamefillfactor',[pixelSize pixelSize]);   
pixel = pixelSet(pixel,'conversiongain', conversiongain);        
pixel = pixelSet(pixel,'voltageswing',voltageSwing);                                             
pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;               
pixel = pixelSet(pixel,'readnoisevolts',readnoise);  
pixel = pixelSet(pixel,'wave',wave);
pixel = pixelSet(pixel,'pixelspectralqe',qe); 


% Set these sensor properties
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'autoExposure','off'); 
sensor = sensorSet(sensor,'rows',rows);
sensor = sensorSet(sensor,'cols',cols);
sensor = sensorSet(sensor,'dsnulevel',dsnu);  
sensor = sensorSet(sensor,'prnulevel',prnu); 
sensor = sensorSet(sensor,'analogGain',analogGain);     
sensor = sensorSet(sensor,'analogOffset',analogOffset);
sensor = sensorSet(sensor,'quantizationmethod',quantizationMethod);

% Stuff the pixel back into the sensor structure
sensor = sensorSet(sensor,'pixel',pixel);
sensor = pixelCenterFillPD(sensor,fillfactor);

% Then we load the calibration data and attach them to the sensor structure
sensor = sensorSet(sensor,'filterspectra',ones(nWaves,1));
sensor = sensorSet(sensor,'infraredfilter',ones(nWaves,1));


%% Optics

oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'model','DiffractionLimited');
optics = opticsSet(optics,'off axis method','Skip');
optics = opticsSet(optics,'focallength',focalLength);
optics = opticsSet(optics,'fnumber',fNumber);

end

