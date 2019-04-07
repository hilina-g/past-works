%% 14 July 2016 Hilina Gudeta, Castellanos Aguirre
%% Pocessing image from a video
%% Plot the waves of every frame 
%% Surface tension project 
% This function can plot smoothed waves of a video, makes an animation from
% the plots, gets the wavelength of each frame, and call another funtion which fits the movement of the speaker
% (input excitation of the experiment) with a sinusoidal.
function [mnyy,wavelength,norAM] = plotEachFrame3(NumFr)
% video = VideoReader('square0870_soap_55Hz_110_800fps.avi');
% NumFr = video.NumberOfFrames;
% fps = 800; % velocity of the camara 
% fq = 55; %frequency of the exitation in Hz
mnyy = zeros(NumFr,1);
lpr = 1045/0.05; % length/pixels ratio in this case 50 mm = 605 pixels
maxPvalue = zeros(NumFr,1);
maxIndex = zeros(NumFr,1);
maxPvalue2 = zeros(NumFr,1);
maxIndex2 = zeros(NumFr,1);
maxA = zeros(NumFr,1);
% peakLoc = zeros(NumFr,1);
wavelength = zeros(NumFr,1);
norAM = zeros(NumFr,1);
% fileformatSpec = '%04.0f'; %format numbers with 4 digits / Used on the loop to get an animation from all frames
% mkdir 'output_movie'; %create folder to save the images that were formed
% figure(2);
for im = 1:NumFr
    supo = imread(['video_frames/frame', num2str(im),'.jpg']);% read the desired image gpuArray()
    supcrop1 = imcrop(supo,[225 200 870 200]); % roughly remove the frame of the cell case 
    sup1 = imbinarize(supcrop1,'adaptive','ForegroundPolarity','dark','Sensitivity',0.5);  % Turns sup2 into image into a binary image and define the thershold value to minimize the interclass variance of the black and white 
    sup2 = ~sup1; % inverts the black and white: white wave, black background
    supB = bwconncomp(sup2); % returns the connected components supB found in the binary image sup2
    numPixels = cellfun(@numel,supB.PixelIdxList); % counts the number of pixels  for each supB (group of white pixels)
    [biggest] = max(numPixels); % selects the largest group of white concatenated pixels
    supL = bwareaopen(sup2,biggest); % removes all connected components (objects) that have fewer than "biggest" pixels from the binary image sup2, producing another cleaned binary image, supL.
    supmaj = bwmorph(supL,'majority'); % Sets a pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1s; otherwise, it sets the pixel to 0.
    supA = ~supmaj; % inverts the white wave into a black one and also inverts the background
    supedge = edge(supA,'sobel'); % Get the edge ( upper and lower waves) in a one pixel line
    supcrop2 = imcrop(supedge,[5 1 664 195]); % Crops the image and give a smaller area to work with
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
    cMeters= (1/lpr)*c;        %Change from pixels to milimiters
    meanyyMeters= (1/lpr)*meanyy;%Change from pixels to milimiters
    maxA(im,1)= (max(meanyyMeters)- min(meanyyMeters))*1000;

    
    %% First subplot 
%       Bin image
    subplot(3,2,[2 1])
    %
    plot(cMeters,meanyyMeters,'LineWidth',3,'color','k'); % plots a beatiful ready to process image of each wave
    %
    title(['Max amplitude  :',num2str(maxA(im,1)),'mm']);
    axis([0 0.031 -0.002 0.002]); % fix the frame in the original size
    xlabel('m')
    ylabel('hight (m)') % 'Orientation','horizontal'
    
%% Second subplot 
    % Croped Video
    subplot(3,2,[4 3])
    supcrop3 = imcrop(supcrop1,[215 50 880 200]);
    imshow(supcrop3)

 %% Fifth subplot and sixth
    
    %FFT     
     subplot(3,2,5 )
    m = length(meanyyMeters);% %%Sample length%% in the space domain signal(function) our points which have the sort of sinusiodal form also can be seen as % Window length%
    n = 1000 * 2^nextpow2(m);% New length to a one that fits better the FFT    
    fftmnyyM = fft(meanyyMeters,n);% Discret Fourier Transform (DFT), actually is a Fast Fourier Transform meanyy instead of meanyyMeters
    fs = lpr; % Sampling frequency in meters obtained NSamples(pixel) * lpr
    xrange = (0:n/2-1)*(fs/n); % range of the X axis
    MEANYY_Mag = abs(fftmnyyM);  
    Meanyy_Mag = MEANYY_Mag(1:length(MEANYY_Mag)/2);% gets the absolute value of the Fourier Transform
    [maxPvalue2(im,1) , maxIndex2(im,1)]= max(Meanyy_Mag); % gets the index of the maximum as well as its value
    xmarkers2 =xrange(maxIndex2(im)); % place markers at these x-values
    ymarkers2 = Meanyy_Mag(maxIndex2(im,1));% place markers at these y-values
    
    peakLoc = 1/xrange(maxIndex2(im))*(1000);    
    wd = ceil(peakLoc);  
     if (20 > wd) && (wd > 17)        % if wavelength bigger than 17mm
            
             subplot(3,2,6)
            xmarkers2 =xrange(maxIndex2(im)); % place markers at these x-values
            ymarkers2 = Meanyy_Mag(maxIndex2(im,1));% place markers at these y-values
             
             plot(xrange,Meanyy_Mag,xmarkers2,ymarkers2,'ko') % Plots the likelihood to have some wavelength
            
            title('Fast Fourier Transform')
            axis([0 fs/30 0 0.5]);  % divide (fs/30)
            xlabel('Frequency (1/m)')
            ylabel('Power')
            
            % Wavelength
     
            [maxPvalue2(im,1) , maxIndex2(im,1)]= max(Meanyy_Mag); % gets the index of the maximum as well as its value
            wavelength(im) = 1/xrange(maxIndex2(im))*1000; %takes the value of the wavelength (mm)
            subplot(3,2,5)
           
            xpeak = 0:NumFr-1;
            plot(xpeak,wavelength,'*') % plots the value of the readed wavelengths
            axis([0 NumFr 0 50])
            title(['Wavelength  :',num2str(wavelength(im,1)),'mm']);
            xlabel('Frames')
            ylabel('Wavelength (mm)') 
         else
                          %For smaller amplitudes
             subplot(3,2,6)
             Meanyy_Mag2 = Meanyy_Mag(2000:length(Meanyy_Mag));%  Meanyy_Mag2 = Meanyy_Mag(7500:length(Meanyy_Mag));
             xrange2= xrange(2000:length(Meanyy_Mag));   %xrange2= xrange(7500:length(Meanyy_Mag));
             [maxPvalue(im,1) , maxIndex(im,1)]= max(Meanyy_Mag2);
             xmarkers =xrange2(maxIndex(im)); % place markers at these x-values
             ymarkers = Meanyy_Mag2(maxIndex(im,1));% place markers at these y-values
        %
             plot(xrange2,Meanyy_Mag2,xmarkers,ymarkers,'ko')
             title('Fast Fourier Transform')
             axis([0 fs/30 0 0.5]);  % divide (fs/30)
             xlabel('Frequency (1/m)')
             ylabel('Power')
                         % To compare
%                         subplot(3,2,3)
                        [maxPvalue2(im,1) , maxIndex2(im,1)]= max(Meanyy_Mag); % gets the index of the maximum as well as its value
                        wavelength(im) = 1/xrange(maxIndex2(im))*1000; %takes the value of the wavelength
                        xmarkers2 =xrange(maxIndex2(im)); % place markers at these x-values
                        ymarkers2 = Meanyy_Mag(maxIndex2(im,1));% place markers at these y-
%                          plot(xrange,Meanyy_Mag,xmarkers2,ymarkers2,'ko') % Plots the likelihood to have some wavelength
%                          axis([0 fs/30 0 5000]);
%                          xlabel('Frequency (1/m)')
%                          ylabel('Power')
% %          Wavelength  
            
            [maxPvalue(im,1) , maxIndex(im,1)]= max(Meanyy_Mag2); % gets the index of the maximum as well as its value
            wavelength(im) = 1/xrange2(maxIndex(im))*1000; %takes the value of the wavelength
             subplot(3,2,5)
           
            xpeak = 0:NumFr-1;
            plot(xpeak,wavelength,'*') % plots the value of the readed wavelengths
            axis([0 NumFr 0 50])
            title(['Wavelength  :',num2str(wavelength(im,1)),'mm']);
            xlabel('Frames')
            ylabel('Wavelength (mm)') 
     end
   pause(0.001);
   for am = 1:NumFr
   norAM(am) = 1/wavelength(am)*maxA(am);
   end
end
