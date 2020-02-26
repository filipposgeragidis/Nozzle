clear all
clc

mdot = 0.1367;
P = 25e5;
T= 3252;
g = 1.1726;
R = 318.6;


A = ((mdot * sqrt(T)) / P) * (sqrt(R/g)) / (((g+1)/2)^(-(g+1)/(2*(g-1))))

rad = sqrt(A/pi)