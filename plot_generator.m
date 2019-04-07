%% 14th July 2016 Castellanos Aguirre and Hilina Gudeta
%% Pocessing image from a video
%% Plot the waves of every frame 
%% Surface tension project 

% Plot_generator can plot smoothed Faraday waves, the original video from 
% those amazing waves were extracted, the teachings of Jean-Baptiste Joseph Fourier 
% as a likelyhood of a frequency to happen and the dominant wavelength all
% of these for each frame you want to process. 

% You will also recieve three values, one is the wavelength for each frame 
% in an array format, while the second will be "mnyy" value, which translates
% the Faraday wave (FW) into dots in a 2D space for MATLAB, this vector will be used in
% the "fitting_funtion_for_the_plot_generator". The third value is the normalization of a
% length, which will allow you to have a third paramater in the 2D plot.

% This actual function needs a number of frames "NumFr" as a reference to work with

function [mnyy,wavelength,norAM] = plot_generator(NumFr,Fs)  
mnyy = zeros(NumFr,1); % Prealocation of the FW
%  Fs % pixels/length   e.g. pixels/ 1130-99% / 1045-70%
maxPvalue = zeros(NumFr,1);  % Prealocation of a vector
maxIndex = zeros(NumFr,1);   % Prealocation of a vector
maxA = zeros(NumFr,1);       % Prealocation of a vector
wavelength = zeros(NumFr,1); % Prealocation of a vector
norAM = zeros(NumFr,1);      % Prealocation of a vector
%% The next two commands and two more at the bottom will allow you to record images of the frames you want.
% You can comment them in order to not use them. At the bottom of the
% script is the second part of this record.
% fileformatSpec = '%04.0f'; %format numbers with 4 digits / Used on the loop to get an animation from all frames
% mkdir 'output_movie'; %create folder to save the images that were formed

%%
figure(1);   % Within this figure you will see ploted al the plots you customize 
for im = 1:NumFr   % For loop which processes all the frames, one by one
    supo = imread(['video_frames/frame', num2str(im),'.jpg']);% reads the frame
    supcrop1 = imcrop(supo,[200 200 850 200]); % roughly removes the frame of the cell case It depends on the video e.g. % 70-[125 200 870 200]/ eth 95%[200 200 850 200] / soapy water [350 250 660 200]
    sup1 = imbinarize(supcrop1,'adaptive','ForegroundPolarity','dark','Sensitivity',0.5);  % Turns sup2 into image into a binary image and define the thershold value to minimize the interclass variance of the black and white 
    sup2 = ~sup1; % inverts the black and white: white wave, black background
    supB = bwconncomp(sup2); % returns the connected components supB found in the binary image sup2
    numPixels = cellfun(@numel,supB.PixelIdxList); % counts the number of pixels  for each supB (group of white pixels)
    [biggest] = max(numPixels); % selects the largest group of white concatenated pixels
    supL = bwareaopen(sup2,biggest); % removes all connected components (objects) that have fewer than "biggest" pixels from the binary image sup2, producing another cleaned binary image, supL.
    supmaj = bwmorph(supL,'majority'); % Sets a pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1s; otherwise, it sets the pixel to 0.
    supA = ~supmaj; % inverts the white wave into a black one and also inverts the background
    supedge = edge(supA,'sobel'); % Get the edge ( upper and lower waves) in a one pixel line
    supcrop2 = imcrop(supedge,[1 1 849 200]); % Crops the image and give a smaller area to work with   /Soapy water [1 10 650 110]
    %%dividing the image
    [r, c] = find(supcrop2==1); %     cr = [c, r];
    [p, f] = find(supcrop2==1); %     fp = [f, p];
    %% Higher curve
    repeatedvaluesC = find(diff(c)==0); % identifies repeated values in the ordinate axis
    c (repeatedvaluesC) = []; % cleans repeated values in the abscissa axis
    r (repeatedvaluesC) = []; % cleans repeated values in the ordenate axis
    %% Lower curve
    p (repeatedvaluesC+1) = [];% clean repeated values in the ordenate axis
    %% Mean curve
    yyp = smooth(p,'moving'); % smoothes the p vector
    yyr = smooth(r,'moving'); % smoothes the r vector
    q = zeros(length(c),1);
    for n = 1:length(c)
        q(n) = (yyp(n,1)+yyr(n,1))/2; % gets the mean value for the ordenate
    end
    yy = smooth(q,'moving'); % smoothes the q vector 
    mnyy(im,1) = mean(yy);  % obtains the mean value of the smoothed vector yy
    meanyy = yy-mean(yy);  % deduct the mean value of the vector yy, which is the movement of the speaker
    cMeters= (1/Fs)*c;      %Change from pixels to milimiters
    meanyyMeters= (1/Fs)*meanyy; %Change from pixels to milimiters
    maxA(im,1)= (max(meanyyMeters)- min(meanyyMeters))*1000; % Gets the Max amplitude of each frame

 %% First subplot   
%    % Bin image
     subplot(3,2,[1 2]) % first row of the plot
%     %
      plot(cMeters,meanyyMeters,'LineWidth',3,'color','k'); % Plots the binarized wave
%     % 
    title(['Max amplitude  :',num2str(maxA(im,1)* 1000),'mm']); % Gives the max amplitude per frame
    xlabel('m')
    ylabel('hight (m)') % 'Orientation','horizontal'
    axis([0 0.029 -0.002 0.002]);%% values for the 605 pix video [0 0.041 -0.002 0.002]
%    
  % %% Second subplot 
% %    % Croped Video
    subplot(3,2,[3 4]) % Second row of the plot
    supcrop3 = imcrop(supcrop1,[1 1 849 200]);% Plots the wave but without much processing, just a crop
    imshow(supcrop3)

 %% Fifth subplot and sixth
    
    %FFT Preparing the Fast Fourier Transform  
    
    m = length(meanyyMeters);% %%Sample length%% in the space domain signal(function) our points which have the sort of sinusiodal form also can be seen as % Window length%
    n = 1000 * 2^nextpow2(m);% New length to a one that fits better the FFT    
    fftmnyyM = fft(meanyy,n);% Discret Fourier Transform (DFT), actually is a Fast Fourier Transform
    fs = Fs; % Sampling frequency in meters obtained NSamples(pixel) * lpr
    xrange = (0:n/2-1)*(fs/n); % range of the X axis
    MEANYY_Mag = abs(fftmnyyM);  % Gets the magnitude of the complex number given by the 
    Meanyy_Mag = MEANYY_Mag(1:length(MEANYY_Mag)/2);% gets the absolute value of the Fourier Transform
    [maxPvalue(im,1) , maxIndex(im,1)]= max(Meanyy_Mag); % gets the index of the maximum as well as its value
        % Fast Fourier Transform     
             subplot(3,2,5) % third row and first column of the plot
             xmarkers =xrange(maxIndex(im)); % place markers at these x-values
             ymarkers = Meanyy_Mag(maxIndex(im,1));% place markers at these y-values
             plot(xrange,Meanyy_Mag,xmarkers,ymarkers,'ko')
             title('Fast Fourier Transform')
             axis([0 fs/30 0 3100]);  % divide (fs/30)
             xlabel('Frequency (1/m)')
             ylabel('Power')
         % Wavelength  
            wavelength(im) = 1/xrange(maxIndex(im))*1000; %takes the value of the wavelength
            subplot(3,2,6)
            xpeak = 0:NumFr-1;  
            plot(xpeak,wavelength,'*') % plots the value of the readed wavelengths
            axis([0 NumFr 0 10])
            title(['Wavelength  :',num2str(wavelength(im,1)),'mm']);
            xlabel('Frames')
            ylabel('Wavelength (mm)') 
                       
     pause(0.001); % Gives a fraction of time to plot each frame
%% This two commands are the second part of the frame recording 
%    F = getframe(figure(1)); %screenshot of the active figure
%    imwrite(F.cdata,['output_movie/frame_',num2str(im,fileformatSpec),'.tif']);
   
     end
   
   for am = 1:NumFr % Just choose one of the following lines or write a new one
%        norAM(am) = 1/max(maxA)*maxA(am);% Normalization given by the Amplitude of eachframe and the MaxAmplitude in the video
   norAM(am) = 1/wavelength(am)*maxA(am); %% Normalization driven by the Wavelength of each video and the correspondent Amplitude
   end
end
