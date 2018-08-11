%--"A model for scattered pressure waveforms from large amplitude oscillations of coated bubbles"--

%--Sijl, Jeroen, Emmanuel Gaud, Peter J. A. Frinking, Marcel Arditi, Nico De Jong, Detlef Lohse, and Michel Versluis. 
%"Acoustic Characterization of Single Ultrasound Contrast Agent Microbubbles." The Journal of the Acoustical Society of America 124.6 (2008): 4091. Print.



function Waveform = ScatteredPressureWaveform(x1,x2,x3)

%---INPUT PARAMETERS---
% t:        time, (unit: s)
% x:        [x(0),x(1),x(3)] radius, velocity and acceleration of bubble wall
% Pa:       acoustic pressure amplitude, (unit: Pa)
% f:        acoutic driving frequency, (unit: Hz)
% r:        distance from bubble to receiving transducer, (unit: m) **NEW**
%---Bubble PARAMETERS---
% R0:       equilibrium bubble radius, (unit: m)
% P0:       hydrostatic pressure (ambient pressure), (unit: Pa), (1.01e5 Pa)
% GAMMA:    polytropic exponent of the gas in the bubble, (1.07 for C3F8)
% MU:       shear (dynamic) liquid viscosity, (unit: Pa), (0.001 Pa.s for water)
% SIGMA:    surface tension, (unit: N/m), (0.072 N/m for water/air) 
% RHO:      density of the liquid, (unit: kg/m3), (998 kg/m3 for water)
% c:        speed of the sound in the liquid, (unit: m/s), (1500 m/s for water)

% CHI:      elasticity parameter of the shell, (unit: N/s)
% KAPPA_S:  shell viscosity, (unit: kg/s)


global RHO r



Waveform = ((RHO*(x1^2)*x3+(2*x1*(x2^2)))/r);


end

