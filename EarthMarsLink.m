clc
clear all
format long

%Global Variables%
global NUM
global BANDS
NUM = 1000;
BANDS = 4;

[dist, t] = distanceCalc();
Rb = linkBudget(dist);
%Plot Things
%Convert time from seconds to Earth Years%
t = t./(365*24*3600);
figure(1);
subplot(3,1,1)
plot(t,dist/10^3,'r');
subplot(3,1,2)
hold;
plot(t,Rb(1,:),'g');
plot(t,Rb(2,:),'b');
subplot(3,1,3)
hold;
plot(t,Rb(3,:),'g');
plot(t,Rb(4,:),'b');

function [dist, t] = distanceCalc()

%CONSTANTS%
%number of iterations (unitless)
global NUM

%GIVENS%
%Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/%
%Start time is 2007 Jan 15, 0030 UT, Earth Time of Periapsis
%Earth is 1, Mars is 2
%Perihelion (m)
r_p = [147.10*10^9, 206.62*10^9];
%Aphelion (m)
r_a = [152.10*10^9, 249.23*10^9];
%Standard Gravitational Parameter of Sun (m^3/s^2)
mu = 1.327124400*10^20;
%Time Range (seconds)
%Counting up from 2007 Jan 15 at 0030 UT.
t = linspace(13*365*24*3600,18*365*24*3600,NUM);

%ORBIT CHARACTERISTICS - from Curtis, Table 3.1%
%Eccentricity (unitless)
ecc = (r_a-r_p)./(r_a+r_p);
%semimajor & semiminor axis (m)
a = (r_p+r_a)./2;
b = a.*sqrt(1-ecc.^2);
%Angular Momentum (m^2/s)
h = sqrt(r_p.*(1+ecc)*mu);
%period (seconds)
T = (((2*pi)/(mu^2)).*(h./(sqrt(1-ecc.^2))).^3);

%VARIABLES%
%True Anomaly (radians)
theta_base = linspace(0,2*pi,NUM);

%CALCULATIONS%
%180deg offset so it starts at perihelion when t = 0.
t_new(1,:) = T(1)/2+t;
%Martian perihelion occurs on 01 Jun 2007 at 0720
%This is 137 days, 6 hours, and 50 minutes after Earth perihelion
t_mars = 137*24*3600 + 6*3600 + 50*60;
t_new(2,:) = T(2)/2+t-t_mars;
%Mean Anomaly-eq3.12 (rad)
for i = 1:2
    M_e(i,:) = 2*pi*t_new(i,:)/T(i);
end
%Solve Keplers Equation%
%Initial guess for E
E = zeros(2,NUM);
for i = 1:2
    if(M_e(i) > pi)
        E(i,:) = M_e(i,:)-ecc(i)/2;
    else
        E(i,:) = M_e(i,:)+ecc(i)/2;
    end
end

% Set the acceptable error
error = 10^(-10);
diff = [1,1];

%iterate until acceptable error is reached
for i = 1:2
    while (diff(i)>error)
        diff(i) = ((E(i)-ecc(i)*sin(E(i))-M_e(i))./(1-ecc(i)*cos(E(i))));
        E(i) = E(i) - diff(i);
    end
end

%Argument of Periapsis + RAAN (rad)
alpha = [102.9*pi/180,(286.5+49.6)*pi/180];
%find the true anomaly
theta = 2*atan((sqrt((1+ecc)/(1-ecc))*tan(E/2)));
%Offset from Solar Center
for i = 1:2
    C_x(i) = (a(i) - r_p(i))*cos(alpha(i));
    C_y(i) = (a(i) - r_p(i))*sin(alpha(i));
end

%Determine complete orbital ellipse
for j = 1:2
    for i = 1:NUM
        x(j,i) = a(j)*cos(alpha(j))*cos(theta_base(i))-b(j)*sin(alpha(j))*sin(theta_base(i))+C_x(j);
        y(j,i) = a(j)*sin(alpha(j))*cos(theta_base(i))+b(j)*cos(alpha(j))*sin(theta_base(i))+C_y(j);
    end
end

%Earth's position at given time
x_E = a(1)*cos(alpha(1))*cos(theta(1,:))-b(1)*sin(alpha(1))*sin(theta(1,:))+C_x(1);
y_E = a(1)*sin(alpha(1))*cos(theta(1,:))+b(1)*cos(alpha(1))*sin(theta(1,:))+C_y(1);
%Mars' position at given time
x_M = a(2)*cos(alpha(2))*cos(theta(2,:))-b(2)*sin(alpha(2))*sin(theta(2,:))+C_x(2);
y_M = a(2)*sin(alpha(2))*cos(theta(2,:))+b(2)*cos(alpha(2))*sin(theta(2,:))+C_y(2);
%L1 Mars Position
%L1 Earth Position

%Distance between Earth and Mars (m)
dist = sqrt((x_E-x_M).^2 + (y_E-y_M).^2);

%Plot Results
%f = figure;
%plot(t,dist,'b');

end

function Rb = linkBudget(dist)
%INFO%
%Calculates max theoretical data rate at some distance d%
%1 is X, 2 is Ka%

%CONSTANTS%
%Number of bands used (unitless)
global BANDS
%Number of iterations on distance
global NUM
%Speed of Light (m/s)
c = 3*10^8;
%Boltsmann's constant (dB(W/k/Hz)
k = -228.5991672;

%GIVENS%
%Frequency (MHz)
f = [8400,32000,1*10^5, 2*10^5];
%Wavelength (m)
lambda = c./(f.*10^6);
%Power (Watts)
P = [5000, 5000,5000,5000];
%Power (dBW)
P_dB = 10*log10(P);
%Transmitter Line Loss (dB, unitless)
LL = 1;
%Transmitter Antenna Diameter (m)
d_Tx = [8.5,8.5,8.5,8.5];
%Transmitter and Receiver Antenna Efficiency (unitless)
eff = 0.65;
%Transmitter Antenna Gain (dB)
Gt = 10*log10(((pi.*d_Tx./lambda).^2).*eff);
%Pointing Loss (dB, unitless)
PL = 1;
%Free Space Loss (dB, unitless)
for i = 1:NUM
    for j = 1:BANDS
        %note, distance needs to be converted to Km%
        FSL(j,i) = 32.4+20*log10(dist(i)/(1*10^3))+20*log10(f(j));
    end
end
%Receiver Antenna Diameter (m)
d_Rx = [8.5,8.5,8.5,8.5];
%Receiver Antenna Gain (dB)
Gr = 10*log10(((pi.*d_Rx./lambda).^2).*eff);
%System Noise Temperature (K) - ASSUMPTION
T = [50,80,100,120];
%Sytem Noise Temp (dBK)
T_dB = k+10*log10(T);
%Required Eb/N0 (dB0
EbN0_req = 1;
% Eb/N0 Margin (dB)
EbN0_margin = 3;

%Calculate Data Rate
PtN0 = zeros(BANDS,NUM);
EbN0 = zeros(BANDS,NUM);
Rb = zeros(BANDS,NUM);
for i = 1:NUM
    for j = 1:BANDS
        %Recieved Pt/N0 (dB-Hz)
        PtN0(j,i) = P_dB(j) - LL + Gt(j) - PL - FSL(j,i) + Gr(j) - T_dB(j);
        %Available Eb/N0 (dB)
        EbN0(j,i) = PtN0(j,i) - EbN0_req - EbN0_margin;
        %Maximum Theoretical Data Rate (Mb/s)
        Rb(j,i) = (10^(EbN0(j,i)/10))/10^6;
    end
end

end

