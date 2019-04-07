%% Hilina Gudeta, Castellanos Aguirre June 27th 2016
%% Finding Wavelength
%% Surface tension 

%extracting wave array from frame
supo = imread('frame150.jpg');% read the desired image gpuArray()
sup = supo(360:450,330:880,:);
sup1 = imbinarize(sup,'adaptive','ForegroundPolarity','dark','Sensitivity',0.5);  % Turns sup2 into image into a binary image and define the thershold value to minimize the interclass variance of the black and white 
supthin1 = bwmorph(sup1,'thin');
supthin2 = bwmorph(supthin1,'thin');
supthin3 = bwmorph(supthin2,'thin');
supthin4 = bwmorph(supthin3,'thin');
supthin5 = bwmorph(supthin4,'thin');
supthin6 = bwmorph(supthin5,'thin');
supthin7 = bwmorph(supthin6,'thin');
supmaj = bwmorph(supthin7,'majority');
supSkel = bwmorph(supmaj,'skel',Inf); 
supedge = edge(supmaj,'sobel');
supedgeP1 = supedge(9:84,9:543,:);
%% Higher curve
supedgeP2 = supedge(9:34,13:543,:);
[r, c] = find(supedgeP2==1);
cr = [c, r];
repeatedvaluesC = find(diff(c)==0);
c (repeatedvaluesC) = [];
r (repeatedvaluesC) = [];
% plot(c,r)
% grid on
% hold on
%% Lower curve
supedgeP3 = supedge(34:60,13:543,:);
[o, d] = find(supedgeP3==1);
do = [d, o];
repeatedvaluesD = find(diff(d)==0);
d (repeatedvaluesD) = [];
o (repeatedvaluesD) = [];
% plot(d,o)
% grid on
% hold on
%% Mean curve
n=1;
for n = 1:length(c)
    q(n) = (o(n,1)+r(n,1))/2;
end
plot(d,q)
title('Wavelength Array');
xlabel('Distance in pixels'); 
ylabel('Amplitude');
grid on
hold on
[maxi,I] = max(q);
[mini,i] = min(q);
amplitude = ((maxi-mini)/2)*(50/605);

figure(1);

% calibFactor = 50/605; %50mm is the measured width of the Hele-Shaw cell using and 608 is the measured pixel equivalent found using imtool(cropped150)
% distance1WLInMm = distance1WLInPixels * calibFactor;

%% fft power series plot to analyse most occuring wavelengths in frame 150


Fs = 605/(50*10^-3);   % Sampling frequency: datapoints or pixels/meter
T = 1/Fs; % Sampling period
sigLen = length (q);
L = sigLen;       % Length of signal, or in other words, length(number) of pixels horizontally
%t = (0:L-1)*T;        % Time vector

X = q; % X is the wave equation or array
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L; % the 0:L/2 part determines the domain boundary or length on the x axis, which sometimes be unnecessarily long
figure(2);
plot(f,P1);
title('Most Dominant Wavelength Analysis');
xlabel('1/wavelength in m^-1'); % this is f. the array shown on top of the power series graph has space pixel in x and amplitude in y directions. The power series graph has space in pixels on the and wavelength on the y directions. In this use of fft, space replaces time, but it is still possible to calculate frequency by doing 1/space in pixels instead of 1/time.
ylabel('Power P1(1/wavelength in m^-1)');
[~, locs] = findpeaks(P1, 'Npeaks',3,'SortStr','descend'); % the highest value returned will be taken as the most accurate wavelength reading
inverseWvInM= f(locs(1));
wvInMm = (1/inverseWvInM)*1000;
fprintf('wavelength is %d mm\n', wvInMm)
