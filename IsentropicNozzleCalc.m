%Nozzle Flow Solver
%This script solves 1D isentropic CD nozzle flow for the Spaceport America
%static test-rocket
%Written by Filippos Geragidis on 1/12/18
%Housekeeping
clear;
clc;

%Defining Variables
%Combustion Chamber Exit Temperature (K)
T0=3252;
%Combustion Chamber Exit Pressure (Pa)
p0=4000000;
%Exit Pressure (pe=patm @4200ft of altitude) (Pa)
pe=87000;
%Cp J/K/kg
Cp0=3547;
%Isentropic parameter
gamma=1.1726;
%Maximum Throttle Mass Flow Rate (kg/s)
mdotmax0=0.122;
%Minimum Throttle Mass Flow Rate (kg/s)
mdotmin0=0.0897;
%Combustion Chamber Exit Gas Constant calculation:
%Cv
Cv0=Cp0/gamma;
%R0
R0=Cp0-Cv0;
%Specific Impulse (s)
Isp=250;
%Maximal and Minimal Thrust Settings:
Fmax=300;
Fmin=200;
%Launch Graviational Acceleration (N/kg)
g0=9.80655;

%Mass Flow Rate Determination
mdotmax0=Fmax/(Isp*g0);
mdotmin0=Fmin/(Isp*g0);

%Combustion Chamber Exit Density calculation:
rho0=p0/(R0*T0);

%Throat Quantities:
%Throat Pressure
pt=p0*(2/(gamma+1))^(gamma/(gamma-1));
Tt=T0*(2/(gamma+1));
rhot=rho0*(2/(gamma+1))^(gamma/(gamma-1));

%Exit Mach Number
Me=(((pe/pt)^((1-gamma)/gamma)-1)*(2/(gamma-1)))^0.5;

%Outlet-Throat Area Ratio (epsilon)
epsilon=(1/Me)*((2+(gamma-1)*Me^2)/(gamma+1))^((gamma+1)/(2*(gamma-1)));

%Throat Area
At=mdotmin0/(p0*((gamma/(R0*T0))*((gamma+1)/2)^((gamma+1)/(1-gamma)))^0.5);
dt=(4*At/pi)^0.5;

%Exit Area
Ae=epsilon*At;
de=(4*Ae/pi)^0.5;




