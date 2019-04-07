%% 21st july 2016 Hilina Gudeta, Castellanos Aguirre
%% Pocessing image from a video
%% Surface tension project
%% 
% Hello dear reader,
 
% this script processes a bunch of videos, in order to graph the behavior of
% the wavelength of a fluid through the variation of the frequency. Feel
% free to use it and to improvit please.

%   It will read video after video in a for loop and process them one by
%   one, ploting a series of wave length values for a given frequency. 
%   Each video will be a record of a given frequency.

%  You have to define the velocity with which the video was recorded   fps
%  You have to write down the excitation frequncy of each video 
%  You have to give the pixel/meters frequency  Fs

%  Ass far as this code goes, it will work better for a same angle and
%  distance for al the videos. Take care about that and carefull with the
%  "fps", "fq" and "Fs" that needed to be fixed. 

files = dir('/Users/jesuscastellanos/Documents/MATLAB/Processing images/Atomate model/square004_scotch_eth95_800fps/square004_scotch_eth95_*.avi');
% makes an array with all the square004_scotch_eth95_*.avi files. Just
% follow the format and give the location of the files
file = files'; % Transpose the files vector

for listOfFiles=1:12          % Reads each video through the list within a for loop
    theVideo =(file(listOfFiles).name)  % Select the video from the array "files"
    video = VideoReader(theVideo); % turns the .avi file into a MATLAB recognizable file type to work with 
    NumFr = video.NumberOfFrames; % gives you the number of frames that has the selcted video
    fps = 800; % Frames per second: velocity of the camera 
    fq= [100,110,110,112,45,50,55,60,65,70,75,80,85,85,90]; %frequency of the exitation in Hz of each video on the list. the number with 1(one) as first digit will be the first in the row
    % for ethanol_70 and the listOfFiles=38 [100, 105, 110, 115, 32, 33, 34, 35, 37, 38, 39, 40, 40, 41, 42, 42, 43, 44, 45, 48, 48, 50, 52, 53, 54, 55, 56, 58, 60, 65, 70, 75, 76, 78, 80, 85, 90, 95];    
     Fs = 1119/(50*10^-3); % Size ratio [pixels/m]  1130-99% / 1045-70%
    % 50 mm is the measured width of the Hele-Shaw cell using and 1130 pixels is the measured pixel 
    % equivalent found using imtool(frameX). 
     mkdir 'video_frames';  % create folder to save the frames that were taken from the video
    
    %% Extract frames from video
    for img = 1:NumFr    % This for loop sweeps through all frames in order to name them and give them file format
        FrameN = strcat('frame', num2str(img),'.jpg'); % Order the frames in a horizontal vector and gives them the file format jpg
               n = read(video,img); %get the values of each frame NumFr
         imwrite(n,['video_frames/frame',num2str(img),'.jpg']); % Gives the name to each frame
    end
    
    %% Plot every wave of the video
    
    [mnyy,wavelength,norAM] = plot_generator(NumFr,Fs); 
    
    %% Getting the mean value of the wavelength
    %  [mwlD,wavelengthOZ] = mean_wave_length(gwavelength);
   
    %% Graph the readings

    figure(2);
    X = ones(size(wavelength))*fq(listOfFiles);
    freq = X(:,1)';
    colormap jet
    scatter(freq,wavelength, [],norAM,'*')
    hold on
    h=colorbar;
   
    xlabel('Excitation Frequency (Hz)')
    ylabel('Wavelength (mm)')
    axis([30 120 0 15])
    
    rmdir ('video_frames','s'); % delete the folder in which were saved all the frames that were taken from the video
end
