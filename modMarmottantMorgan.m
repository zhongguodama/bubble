% --"A model for large amplitude oscillations of coated bubbles accounting for buckling and rupture"--
% --Philippe Marmottant, et al., J. Acoust. Soc. Am. 118(6), 2005--

function rprime = modMarmottantMorgan(t,x) 
%---INPUT PARAMETERS---
% t:        time, (unit: s)
% x:        [x(0),x(1)]_intial,intial radius and velocity of bubble
% Pa:       acoustic pressure amplitude, (unit: Pa)
% f:        acoutic driving frequency, (unit: Hz)
%---Bubble PARAMETERS---
% R0:       equilibrium bubble radius, (unit: m)
% P0:       hydrostatic pressure (ambient pressure), (unit: Pa), (1.01e5 Pa)
% GAMMA:    polytropic exponent of the gas in the bubble, (1.07 for C3F8)
% MU:       shear (dynamic) liquid viscosity, (unit: Pa), (0.001 Pa.s for water)
% SIGMA:    surface tension, (unit: N/m), (0.072 N/m for water/air) 
% RHO:      density of the liquid, (unit: kg/m3), (998 kg/m3 for water)
% c:        speed of the sound in the liquid, (unit: m/s), (1500 m/s for water)

% CHI:      elasticity parameter of the shell, (unit: N/s) - NOTE: The CHI term was taken from Microbubble spectroscopy of ultrasound contrast agents paper 
% KAPPA_S:  shell viscosity, (unit: kg/s)



global CHIX R0 P0 GAMMA MU RHO c CHI KAPPA_S Sigma_initial Rbreakup Rbuckle kexi0 r_rupture Pa f t_end tspan_excitation window_length t_interval hannwin

%% nonlinear shell behavior of coated microbubbles

% surface tension

if x(1)<=Rbuckle
    SIGMA=0;
elseif  (x(1)>Rbuckle&&x(1)<=r_rupture)
        SIGMA=CHI*((x(1)/Rbuckle)^2-1);
    else
        SIGMA=kexi0;
end


SIGMA_R0 = Sigma_initial; 

%% driving pulse, sine burst+hanning win

if  (0<=t&&t<=t_end)
        P_drive = Pa*sin(2*pi*f*t).*hannwin(floor(t/t_interval+1));
   else
       P_drive = 0;
end

%P_drive = [Pa*sin(2*pi*f*tspan_excitation).*hann(length(tspan_excitation))',zeros(1,length(t)-length(tspan_excitation))];

rprime = zeros(2,1);

rprime(1) = x(2);       % velocity
rprime(2) = 1/(RHO*x(1)) * (  (P0+2*SIGMA_R0/R0+2*CHIX/R0)*(x(1)/R0)^(-3*GAMMA)*(1-3*GAMMA/c*x(2)) - P0 - 2*SIGMA/x(1)...
              - 2*CHIX/x(1)*(R0/x(1))^2 - 4*MU*x(2)/x(1) - P_drive  - 4*KAPPA_S*x(2)/x(1)^2  )-3/2*(x(2))^2/x(1);  % acceleration
          
% ode45 calculaltion gives back radius (x(:,1))and velocity (x(:,2))