% M-file to calculate the linear extinction cross-section and attenuation
% of a given microbubble population
%
% Based on the paper "uniform scattering and attenuation..."
% fit 2kpa attenuation

clear all 
%close all

tic
%% load size distribution and exp data
load size20171214tip730ul
load 2017-12-14dspc2_3ul
% load size20171214dspc2_3ul.mat % size data
% load 2017-12-14dspc2_3ul.mat     % fexcitation-attenuation data

plott = 1;% set this to 0 when fitting

%% driving parameters
frange = [1:0.1:2]*1e6;% frequency range of interest
Pressure = 40e3; % driving pressure Pa
                    
%% bubble parameters
r0 = 0.5*D(:,1)/1e6;% Microbubble radii of interest (rmin to rmax of sizing data), um
%r0 = linspace(2e-6,3e-6,100);% Microbubble radii of interest (rmin to rmax of sizing data)
% plot mb size distribution
figure;
plot(D(:,1),D(:,2))

elas_mod_range = 3.3%:0.01:0.9%:0.1:1;% Elastic compression modulus (Marmottant model), N/m Chi
kaps_range =  1e-8%[1:1:25]*1e-9;% Interfacial dilational viscosity
sig_mar_range = 0.002%[0:0.001:0.02];% initial Surface tension water/lipid/PFC interface (coated bubble), N/m
nratio = 8% adjustment of mbs concentration

%% holder parameters
volume = 50*1e-6; % holder volume m^3 5e-6 50ml
d = 0.5e-2; % holder path length m

%% Physical Parameters (SI units)

kap = 1.07;% Octofluoropropane polytropic exponent
mu = 0.001;% Shear liquid viscosity, Pa*s
rho = 998;% Density of liquid, kg/m^3
P0 = 1e5;% Atmospheric pressure, Pa
sig = 0.072;% Surface tension, water/PFC gas interface (free bubble), N/m
c = 1498;% Sound speed in liquid, m/s


for CHIii = 1:length(elas_mod_range)
    for kapsii = 1:length(kaps_range)
        for sig_marii =1:length(sig_mar_range)
            elas_mod = elas_mod_range(CHIii);
            kaps = kaps_range(kapsii);
            sig_mar = sig_mar_range(sig_marii);
            
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
                    
                    % Coated bubble (linear Marmottant)
                    cap_omega_mar = w0_mar./wd;
                    lin_ext_cross_mar(iradius,:) = (4*pi*r0(iradius)^2).*((c.*d_total_mar)./(w0_mar*r0(iradius))).* ...
                        ((cap_omega_mar.^2)./(((1-(cap_omega_mar.^2)).^2)+((cap_omega_mar.^2).* ...
                        (d_total_mar.^2))));% Linear extinction cross-section for Marmottant Model
                   
                    % Coated bubble (full Marmottant)
                    for fdriveii = 1:length(frange)
                    [Pd,Ps] = scatterWave(r0(iradius),elas_mod,kaps,sig_mar,frange(fdriveii),Pressure);
                    nonlin_ext_cross_mar(iradius,fdriveii) = 4*pi*r0(iradius)^2*sum(Ps.^2)/sum(Pd.^2);
                    end
                    
                end
            
            
            %% Attenuation
            
            mb_dist = nratio*D(:,2);% Microbubble size distribution (number density)
            %mb_dist = 2.2e4.*gaussmf(r0,[.15e-6 2.5e-6]);% Microbubble size distribution (number density)
            
            for iradius = 1:length(r0)
                % Free bubble
                att(iradius,:) = 10*log10(exp(1)).*lin_ext_cross(iradius,:).*mb_dist(iradius);
                
                % Coated bubble (linear Marmottant)
                att_mar(iradius,:) = 10*log10(exp(1)).*lin_ext_cross_mar(iradius,:).*mb_dist(iradius);
                % Coated bubble (nonlinear Marmottant)
                non_att_mar(iradius,:) = 10*log10(exp(1)).*nonlin_ext_cross_mar(iradius,:).*mb_dist(iradius);

            end
            
            
            %% fitting to exp to calculate error
            
            % volume averaged extinction cross section multiples d (equation 5 in paper)
            % coverted into db/cm
            %plot(frange./1e6,sum(att)/volume*d/(d*1e2),'k'? % free bubble
            atten_c_sim = sum(att_mar)/volume*d/(d*1e2); % coated bubble
            
            atten_c_sim_yy = spline(frange./1e6,atten_c_sim,fexcitation);
            error(CHIii,kapsii,sig_marii) = sum((atten_c_sim_yy(1:end) - (atten_c(1,1:end))-2.8).^2); % 10 means fit from f = 2mhz
            
        end
    end
end

toc

%% find best fit
[v,loc] = min(error(:));
[ii,jj,k] = ind2sub(size(error),loc);
CHI_fit = elas_mod_range(ii)
kaps_fit = kaps_range(jj)
sig0_fit = sig_mar_range(k)


%% Plotting
if plott == 1;
    figure(3);
    plot(frange./1e6,atten_c_sim,'k',fexcitation,atten_c_sim_yy,'^',fexcitation(8:end),atten_c(1,8:end),'o');hold on
    legend('linear simulation','linear simulation-experimentf','experiment');
    ylabel('Attenuation [dB/cm]');xlabel('Frequency [MHz]');
    figformat
    xlim([1 3])
end
%plot(frange./1e6,sum(att)/volume*d/(d*1e2),'k') % free bubble
figure(4);hold on; plot(frange./1e6,sum(non_att_mar)/volume*d/(d*1e2),'k');hold on;
legend(['nonlinear simulation - ',num2str(Pressure/1e3),'kPa']);
ylabel('Attenuation [dB/cm]');xlabel('Frequency [MHz]');
figformat

%% volume fraction of mbs
volume_mbs = sum((4/3*pi*r0.^3).*mb_dist); % m^3
volumeFraction = volume_mbs/(volume)