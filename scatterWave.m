function  [Pd,SPW] = scatterWave(radius,elas_mod,kaps,sig_mar,fd,Pressure)
% calculate bubble scattered waveform

global CHIX R0 P0 GAMMA MU RHO c CHI KAPPA_S Sigma_initial Rbuckle kexi0 r_rupture Pa f t_end tspan_excitation window_length t_interval hannwin


%% Constant Parameters
P0 = 1e5;                    % hydrostatic pressure (ambient pressure), (unit: Pa), (1.01e5 Pa)
GAMMA = 1.07;                   % polytropic exponent of the gas in the bubble, (1.4 for air, 1.07 for perfluorobutane (c4F10) gas)
MU = 0.001;                     % shear (dynamic) liquid viscosity, (unit: Pa), (0.001 Pa for water)
kexi0 = 0.072;%725;%0.0725;    % surface tension, (unit: N/m), (0.072 N/m for water/air)
RHO = 998;                     % density of the liquid, (unit: kg/m3), (998 kg/m3 for water)
c = 1498;                       % speed of the sound in the liquid, (unit: m/s), (1500 m/s for water)

%% bubble parameters
R0 = radius;
CHI =  elas_mod;
KAPPA_S = kaps;
Sigma_initial =  sig_mar;
CHIX = 0;

Rbuckle = R0/sqrt( Sigma_initial/CHI+1);                              % Rbuckle
r_rupture = Rbuckle*(sqrt(kexi0/CHI+1));

%% driving parameters
cycle = 10;
win = 0;
f = fd;
Pa = Pressure;

% Numerical integration
x0 = [R0;0];                          % intial conditions of [x(0),x(1)], radius and velocity of bubble wall
t_end = cycle*1/f;                               % integer * cycles. drivewave end time
fs = 50e6;
t_interval = 1/fs;
tspan_excitation = [0:t_interval:t_end];     % tspan of acoustic excitation
window_length = length(tspan_excitation);     % define hanning window length over tspan of acoustic excitation
n = tspan_excitation/t_interval+1;
if win == 1; % hanning wind
    tukey = tukeywin(window_length,opt);
    hannwin = tukey;
    %hannwin = 0.5*(1-cos(2*pi*n/window_length))'; % hann window for modulating driving pulse
else % no window
    hannwin = n./n;
end

tspan = [0:t_interval:5*t_end];   % interval of integration

    [t,x] = ode15s(@modMarmottantMorgan,tspan,x0);


Pd = Pa*sin(2*pi*f*tspan_excitation);

x1 = x(:,1);  % R - t
x2 = x(:,2);  % v - t


%% Calculating Scattered Pressure Waveform
maxrows = length(t);

i = 1;
time = 0;
x3 = zeros(maxrows,1); % acceleration


while i <= length(t)
    
    input_conditions = x(i,:)';
        results =  modMarmottantMorgan(time,input_conditions); % obtaining rprime(1)(velocity) and rprime(2)(accerleration) from modRayleighPlesset
    
    x3(i) = results(2);                                                        % inserting the results from modRayleighPlesset into our output matrix
    i = i + 1;                                                                 % index to replace the correct row in out
    time = time + t_interval;                                                  % incrementing the time by the time step determined
    
end

i = 1;
SPW = zeros(maxrows,1);

while i <= length(t)
    
    SPW(i) = RHO*(x1(i)*x3(i)+2*x2(i)^2);
    i = i + 1;
    
end

end

