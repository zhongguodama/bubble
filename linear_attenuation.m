% M-file to calculate the linear extinction cross-section and attenuation 
% of a given microbubble population
%
% Based on the paper "A newtonian rheological model for the interface of
% microbubble contrast agents", by Chatterjee and Sarkar, Ultrasound in
% Med. & Biol., Vol. 29, No. 12, 2003
%
% Mark Burgess, marktb@bu.edu

clear all%close all

%% Input Parameters

frange = linspace(1e6,4e6,100);% frequency range of interest
r0 = linspace(2e-6,3e-6,100);% Microbubble radii of interest (rmin to rmax of sizing data)

elas_mod = 0.6;% Elastic compression modulus (Marmottant model), N/m Chi
kaps = 1e-9;% Interfacial dilational viscosity

%% Physical Parameters (SI units)

kap = 1.07;% Octofluoropropane polytropic exponent
mu = 0.001;% Shear liquid viscosity, Pa*s
rho = 998;% Density of liquid, kg/m^3
P0 = 1e5;% Atmospheric pressure, Pa
sig = 0.072;% Surface tension, water/PFC gas interface (free bubble), N/m
sig_mar = 0.025;% Surface tension water/lipid/PFC interface (coated bubble), N/m
c = 1498;% Sound speed in liquid, m/s

for ifreq = 1:length(frange)
for iradius = 1:length(r0)
%% Driving, Natural, and Resonant Frequencies

wd = (2*pi).*(frange);% Angular driving frequency
w0 = (1./r0(iradius)).*sqrt(((3*kap*P0)/rho).*(1+((2*sig)./(P0.*r0(iradius))))- ... 
    ((2*sig)./(rho.*r0(iradius))));% Angular Minnaert natural frequency including surface tension (free bubble), pg.183 Acoustic Bubble

% Marmottant
w0_mar = (1./r0(iradius)).*sqrt(((3*kap*P0)/rho).*(1+((2*sig_mar)./(P0.*r0(iradius))))- ... 
    ((2*sig_mar)./(rho.*r0(iradius)))+((4*elas_mod)./(r0(iradius).*rho)));%Angular natural frequency for coated bubble (Marmottant model)

%% Damping

% Damping for free bubble
d_vis = (4*mu)./(rho.*w0.*(r0(iradius).^2));% Damping from liquid (viscosity)
d_shell = 0;% Damping from shell encapsulation (dilational viscosity)
d_rad = ((wd.^2).*r0(iradius))./(w0.*c);% Damping from acoustic radiation (frequency dependent!)

d_total = d_vis + d_shell + d_rad;% Total damping

% Damping for Marmottant Model
d_vis_mar = (4*mu)./(rho.*w0_mar.*(r0(iradius).^2));% Damping from liquid (viscosity)
d_shell_mar = (4*kaps)./(rho.*w0_mar.*(r0(iradius).^3));% Damping from shell encapsulation (dilational viscosity)
d_rad_mar = ((wd.^2).*r0(iradius))./(w0_mar.*c);% Damping from acoustic radiation (frequency dependent!)

d_total_mar = d_vis_mar + d_shell_mar + d_rad_mar;% Total damping

%% Extinction Cross-section

% Free bubble
cap_omega = w0./wd;
lin_ext_cross(iradius,:) = (4*pi.*r0(iradius).^2).*((c.*d_total)./(w0.*r0(iradius))).* ... 
    ((cap_omega.^2)./(((1-(cap_omega.^2)).^2)+((cap_omega.^2).*(d_total.^2))));% Linear extinction cross-section for free bubble

% Coated bubble (Marmottant)
cap_omega_mar = w0_mar./wd;
lin_ext_cross_mar(iradius,:) = (4*pi*r0(iradius)^2).*((c.*d_total_mar)./(w0_mar*r0(iradius))).* ... 
    ((cap_omega_mar.^2)./(((1-(cap_omega_mar.^2)).^2)+((cap_omega_mar.^2).* ... 
    (d_total_mar.^2))));% Linear extinction cross-section for Marmottant Model

end
end

%% Attenuation

mb_dist = 5e6.*gaussmf(r0,[.15e-6 2.5e-6]);% Microbubble size distribution (number density)
for iradius = 1:length(r0)
% Free bubble    
att(iradius,:) = 10*log10(exp(1)).*lin_ext_cross(iradius,:).*mb_dist(iradius);

% Coated bubble (Marmottant)
att_mar(iradius,:) = 10*log10(exp(1)).*lin_ext_cross_mar(iradius,:).*mb_dist(iradius);
end

%% Plotting
figure(1)
plot(frange./1e6,sum(att),'k',frange./1e6,sum(att_mar),'r');hold on
legend('free bubble','coated bubble');
set(gca,'FontSize',16);grid on;
ylabel('Attenuation');xlabel('Frequency [MHz]');
