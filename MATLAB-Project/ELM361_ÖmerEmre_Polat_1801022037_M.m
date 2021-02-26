%Ömer Emre POLAT
%1801022037
%Analog Proje 1

clear all; close all; clc %clear variables close plots and clear console

%///////////SETTING VARIABLES///////////
fsample = 20000; %sample frequency
fm = 50; %highest frequency of m
fc = 250; %highest frequency of c
tm = 1/fm; %period of m
tc = 1/fc; %period of c

%///////////SETTING X AXIS VARIABLES///////////
t = linspace(-1/2,1/2,fsample); %vector of t
n = length(t); %signal length
f = linspace(-n/2, n/2-1, fsample) * (fsample/n); %frequency vector

%///////////SETTING Y AXIS VARIABLES///////////
m = (10*cos((fm/2)*2*pi*t) + 20*cos(fm*2*pi*t)); %vector of m
c = (100*cos(fc*2*pi*t)); %vector of c
y = m .* c; %modulated vector
d = (1*cos(fm*2*pi*t)); %demodulator vector
yd = y .* d; %demodulated vector

%///////////PLOTTING m(t)///////////
plot(t, m); %plot t and m vectors
hold on; %hold on so that plot does not disappear
title("m(t)");
xlabel("t");
ylabel("m(t)");
axis([-tm tm (min(m)-1) (max(m)+1)]); %graph scaling

%///////////CALCULATING M(f)///////////
mfft = fftshift(fft(m));
mfft = abs(mfft/length(m));

%///////////PLOTTING M(f)///////////
figure(); %open new plot figure
plot(f, mfft); %plot f and mfft vectors
hold on; %hold on so that plot does not disappear
title("M(f)");
xlabel("f");
ylabel("M(f)");
axis([(-fm-10) (fm+10) (min(mfft)-1) (max(mfft)+1)]); %graph scaling

%///////////PLOTTING Y(f)///////////
figure(); %open new plot figure
plot(t, y); %plot t and y vectors
hold on; %hold on so that plot does not disappear
title("y(t)");
xlabel("t");
ylabel("y(t)");
axis([-tm tm (min(y)-1) (max(y)+1)]); %graph scaling

%///////////CALCULATING Y(f)///////////
yfft = fftshift(fft(y));
yfft = abs(yfft/length(y));

%///////////PLOTTING Y(f)///////////
figure(); %open new plot figure
plot(f, yfft); %plot f and yfft vectors
hold on; %hold on so that plot does not disappear
title("Y(f)");
xlabel("f");
ylabel("Y(f)");
axis([(-(fm+fc)-50) ((fm+fc)+50) (min(yfft)-50) (max(yfft)+50)]); %graph scaling

%///////////PLOTTING yd(t)///////////
figure(); %open new plot figure
plot(t, yd); %plot t and yd vectors
hold on; %hold on so that plot does not disappear
title("yd(t)");
xlabel("t");
ylabel("yd(t)");
axis([-tm tm (min(yd)-100) (max(yd)+100)]); %graph scaling

%///////////CALCULATING Yd(f)///////////
ydfft = fftshift(fft(yd));
ydfft = abs(ydfft/length(yd));

%///////////PLOTTING Yd(f)///////////
figure(); %open new plot figure
plot(f, ydfft); %plot f and ydfft vectors
hold on; %hold on so that plot does not disappear
title("Yd(f)");
xlabel("f");
ylabel("Yd(f)");
axis([(-(fm+fc)-100) ((fm+fc)+100) (min(ydfft)-50) (max(ydfft)+50)]); %graph scaling

%///////////CALCULATING AGF///////////
lpfilter = @(x) (1.*((-fm) <= x & x <= (fm))); %low pass filter function
agf = lpfilter(f); %filtered values vector

%///////////USING AGF///////////
zfft = ydfft .* agf; %demodulated signal passing through the AGF
z = ifftshift(ifft(zfft)); %inverse fourier of the zfft to z

%///////////PLOTTING z(t)///////////
figure(); %open new plot figure
plot(t, z); %plot t and z vectors
hold on; %hold on so that plot does not disappear
title("z(t)");
xlabel("t");
ylabel("z(t)");
axis([-tm tm (min(z)) (max(z))]); %graph scaling

%///////////PLOTTING Z(f)///////////
figure(); %open new plot figure
plot(f, zfft); %plot f and zfft vectors
hold on; %hold on so that plot does not disappear
title("Z(f)");
xlabel("f");
ylabel("Z(f)");
axis([(-(fm+fc)-100) ((fm+fc)+100) (min(zfft)-50) (max(zfft)+50)]); %graph scaling