%% Hilina Gudeta July 19th 2016
%% determining the surface tension

w=input('what is the frequency of vibration of the liquid (excitation frequency/2 ?) ');
g=9.81;
wv=input('what is the wavelength read in mm?');
k=(2*pi)/wv;
d=input ('what is the density of the liquid in kg/m^3?');

T=((((w)^2)+((g)*(k)))*(d))/((k)^3);
display([num2str(T*1000), ' mN/m'])

