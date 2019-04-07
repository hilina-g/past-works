%% 30 june 2016 Castellanos Aguirre, Hilina Gudeta
%% Pocessing image from a video
%% Surface tension 
%% 
% This script obtains the waves formed in the Helle Shaw cell by means of the
% use of a couple of functions 

files = dir('/Users/hilina_baby/Documents/MATLAB/square004_scotch_eth70*.avi');
%files = dir('/Users/hilina_baby/Documents/MATLAB/square004_scotch_soap_*.avi');
%files = dir('/Users/hilina_baby/Documents/MATLAB/square004_scotch_eth95*.avi');

file = files';

%for listOfFiles=1:38
   % theVideo=(file(listOfFiles).name)

video = VideoReader('square004_scotch_eth70_110Hz_110_800fps.avi'); % reads the video, one has to give the its name
NumFr = video.NumberOfFrames; % Tells you the number of frames that have been taken from the video 
fps = 800; % Frames per second: velocity of the camara 

%fq= [100, 105, 110, 115, 32, 33, 34, 35, 37, 38, 39, 40, 40, 41, 42, 42, 43, 44, 45, 48, 48, 50, 52, 53, 54, 55, 56, 58, 60, 65, 70, 75, 76, 78,80, 85, 90, 95]; % eth70 38 videos excitation frequencies

%fq= [36,40,43,45,48,50,53,55,60,63,65,70,75,80,83,85,90,93]; %soap 19 videos excitation frequencies

%fq= [100,110,110,112,45,50,55,60,65,70,75,80,85,85,90]; %eth90 15 videos excitation frequencies

Fs = 1045/(50*10^-3); % Sampling frequency: datapoints or pixels/m. calibFactor = 50/605; 
% 50mm is the measured width of the Hele-Shaw cell using and 605 is the measured pixel 
% equivalent found using imtool(cropped150). distanceInMm = distanceInPixels * calibFactor.
 
%the measured width of the Hele-Shaw cell using and 605 is the measured pixel equivalent found using imtool(cropped150). distanceInMm = distanceInPixels * calibFactor.
% mkdir 'video_frames'; % create folder to save the frames that were taken from the video
%% Extract frames from video

% for img = 1: NumFr
%     FrameN = strcat('frame', num2str(img),'.jpg'); % Order the frames in a horizontal vector and gives them the file format jpg
%            n = read(video,img); %get the values of each frame NumFr
% %      imwrite(n,FrameN) %write de values of n into a file with the name FrameN + serial number
%      imwrite(n,['video_frames/frame',num2str(img),'.jpg']);
% end
%% Plot every the waves of the video
[mnyy,wavelength,norAM] = plotEachFrame3(NumFr);
%% Fitting the movement of the input excitation to a sinusoidal shape
% [C] = function_fitting(NumFr,fps,fq,mnyy);
%% Getting the mean value of the wavelength
%  [mwlD,wavelengthOZ] = mean_wave_length(gwavelength);

%  figure(2);
% freq = ones(size(wavelength))*fq(listOfFiles);
% colormap jet
% scatter(ones(size(wavelength))*fq(listOfFiles),wavelength, [],norAM,'*')
% hold on
% h=colorbar;
% xlabel('Frequency of Vibration (Hz)')
% ylabel('Wavelength (mm)')
% 
% g = 9806; % [mm/s^2] gravity at sea level
% sigma = 0.00002260; % [N/mm]
% rho = (0.3*1000+0.7*789)./(1000^3); % [kg/mm^3]
% lambda = (2.5:0.01:14);% (1./Fs) 4.784688995215311e-05 multiplying the wavelength vector
% h = 16.5; % [mm] 
% k = 2*pi./(lambda);
% % w = sqrt((k.*tanh(k.*h)).*(g+(sigma./rho)*k.^2));% frequency
% freqThrE = (sqrt(g*2*pi./(lambda)+(sigma./rho)*(2*pi./lambda).^3))./2; %freq is the excitation freq. w in the original equation is half of freq. w is the vibrational frequency. freqThrE=w
% plot(freqThrE,lambda)

%end
% 
% 
