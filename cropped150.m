%% June 22 2016 Hilina Gudeta
%% Pocessing wave image to automate wavelength calculation

I= imread('frame150.jpg');
cropped150= imcrop(I, [290 350 934-290 429-350]);
figure;
imshow(cropped150)
h = imdistline(gca);
api = iptgetapi(h);
pause('on');
pause;
width = api.getDistance();
disp(width)
calibFactor = 50/605; %50mm is the measured width of the Hele-Shaw cell using and 608 is the measured pixel equivalent found using imtool(cropped150)
distance1WLInMm = distance1WLInPixels * calibFactor;

%fast fourier series
Fs = 605/(50*10^-3)   % Sampling frequency: datapoints or pixels/meter
T = 1/Fs;             % Sampling period
L = length (q);       % Length of signal
t = (0:L-1)*T;        % Time vector
     
% h1 = imdistline(gca);
% api = iptgetapi(h1);
% pause('on');
% pause;
% a = api.getDistance(); % a is amplitude
% disp(a)
%   
% h2 = imdistline(gca);
% api = iptgetapi(h2);
% pause('on');
% pause;
% f = api.getDistance(); % f is wave frequency
% disp(f)
X = q; % X is the wave equation or array
Y = fft(q);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);
title('Most Dominant Wavelength Analysis');
xlabel('wavelength (mm)'); % it would be frequency f here, but we have a distance space instead of time.
ylabel('|Power P1(wavelength)|');
