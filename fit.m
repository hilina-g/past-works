%% Hilina Gudeta July 2016

% for eth70
openfig('Amp_gradient .fig');
g=9810; %mm/s^2
d=885.*10^-9; %in kg/mm^3
T=0.0000226; %in N/mm
files = dir('/Users/hilina_baby/Documents/MATLAB/square004_scotch_eth70*.avi');    
freqE= [100, 105, 110, 115, 32, 33, 34, 35, 37, 38, 39, 40, 40, 41, 42, 42, 43, 44, 45, 48, 48, 50, 52, 53, 54, 55, 56, 58, 60, 65, 70, 75, 76, 78,80, 85, 90, 95]; % eth70 38 videos excitation frequencies
file = files';
for listOfFiles=1:38
freqG=freqE(listOfFiles); % in 1/s Freq=excitation frequency, and w=frequency of vibration of the liquid
theVideo=(file(listOfFiles).name)
video = VideoReader(theVideo);
NumFr = video.NumberOfFrames;
%for im=1:NumFr
%wv=wavelength(im);
%k=(2*pi)/wv; % in 1/m
%T=((((freqG)^2)+((g)*((2*pi)/wv)))*(d))/(((2*pi)/wv)^3);
% equation solved for wv becomes:
%wv = 1/3 ((2^(2/3).*pi.* (-2 *d^3* freqG^12* g^3+27 *d^2 *freqG^16* T+3 *sqrt(3)*sqrt(27* d^4 freqG^32* T^2-4 *d^5 *freqG^28 g^3 T))^(1/3))./(d *freqG^6)+(2* 2^(1/3).*pi.* d freqG^2 *g^2)/(-2* d^3 *freqG^12 *g^3+27* d^2 *freqG^16* T+3 *sqrt(3)* sqrt(27 *d^4* freqG^32* T^2-4 *d^5 *freqG^28 *g^3 *T))^(1/3))-(2 .*pi.* g)/(3 freqG^2)) 
%wv= (1/3).*(((-2.*(d^3).*(freqG^12).*(g^3) + 3.*(3^(1/2)).*((27.*(d^4).*(freqG^32).*(T^2)- 4.*(d^5).*(freqG^28).*(g^3).*T)^(1/2)))^(1/3)) + ((2.*(2^(1/3)).*pi.*d.*(freqG^2).*(g^2))/(-2.*(d^3).*(freqG^12).*(g^3) + 3.*(3^(1/2)).*((27.*(d^4).*(freqG^32).*(T^2)- 4.*(d^5).*(freqG^28).*(g^3).*T)^(1/2)))^(1/3))) - (2.*pi.*g/(3.*(freqG^2)))
%wv = 1/3 ((2^(2/3) pi (-2 d^3 freqG^12 g^3+27 d^2 freqG^16 T+3 sqrt(3) sqrt(27 d^4 freqG^32 T^2-4 d^5 freqG^28 g^3 T))^(1/3))/(d freqG^6)+(2 2^(1/3) pi d freqG^2 g^2)/(-2 d^3 freqG^12 g^3+27 d^2 freqG^16 T+3 sqrt(3) sqrt(27 d^4 freqG^32 T^2-4 d^5 freqG^28 g^3 T))^(1/3))-(2 pi g)/(3 freqG^2)
wv = (2.5:0.01:14);% (1./Fs) 4.784688995215311e-05 multiplying the wavelength vector

freqG = (sqrt(g*2*pi./(wv)+(T./d)*(2*pi./wv).^3))./2; %freq is the excitation freq. w in the original equation is half of freq. w is the vibrational frequency. freqThrE=w
plot(freqG, wv)
end

% g = 9806; % [mm/s^2] gravity at sea level
% T = 0.00002260; % [N/mm]
% d = (0.3*1000+0.7*789)./(1000^3); % [kg/mm^3]
% wv = (2.5:0.01:14);% (1./Fs) 4.784688995215311e-05 multiplying the wavelength vector
% 
% 

