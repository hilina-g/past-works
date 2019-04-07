
% show the periodogram of a synthetic data

% synthetic data
fps = 100;                          % Sampling frequency (Hz)
t = 0:1/fps:6;                      %  sec sample

f1=1; % frequency 1
f2=2; % frequency 2

x = 0.8*sin(2*pi*f1*t) ...      % f1 component
+ 1.2*sin(2*pi*f2*(t-2)) ...       % f2 component
+ rand(size(t)); % noise

% data preparation
m = length(x);            % Window length
n = pow2(nextpow2(m)*1);
y = fft(x,n);             % DFT
f = (0:n-1)*(fps/n);      % Frequency range

power = y.*conj(y)/n;     % Power of the DFT

figure;
%plot data
subplot(2,1,1);
plot(t,x,'*'); 
ylabel('value');
xlabel('time (s)');

% periodogram
subplot(2,1,2);
plot(f,power); axis([0 fps/10 0 max(power)]);
xlabel('Frequency (Hz)');
ylabel('Power');
title('{\bf Periodogram}');
