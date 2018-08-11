clear;
%close all;
% ***************************************************
% cd attenuation fold. get attenuation data from water and bubble fold (cd each subfold of water and bubble fold)
% Plot of attenuation specturm for different excitation pressure (different attenuator db)
% Qian
%************************************************************

%% variables
basepath = '/Users/liqian/Google Drive/lab/lab/my papers/lcmb based on narrowband attenuation/narrowband attenuation/';
% fold and file info
date = '2017-12-14';% till 20160425
foldname_1 = 'water'; % water
foldname = 'dspc2_3ul'; % bubble
filename = [date,foldname];
fs = 50e6; % sampling rate in gagescope 
d =3; % cm
plotDetails = 0; % plot voltage-time curve for water and mbs and conresponding power specturm
fig = 101;

% gate_n = {'5000:6249','6250:7499','7500:8749','8750:9999','10000:11249','11250:12499',...
%     '12500:13749','13750:14999','15000:16249','16250:17499'};% gate for each freq
% gate_n = {'5000:6499', '6500:7999','8000:9499','9500:10999','11000:12499','12500:13999','14000:15499'};
subfoldname = {'2kPa','4kPa','6kPa','8kPa','10kPa','20kPa','30kPa','40kPa','50kPa','60kPa','70kPa','80kPa'};
%subfoldname ={'2kPa\','10kPa\','40kPa\','60kPa\','80kPa\'};
%subfoldname = {'10kPa\'}%,'30kPa\','50kPa\','70kPa\','90kPa\'};
file = 3;%4 3-first file in subfoldname(1-3Mhz);4-second file in subfoldname(3-5Mhz);

fp1 = [basepath,date,'/',foldname_1,'/',subfoldname{1},'/'];
dir_result = dir(fp1);
load([fp1,dir_result(file).name]); 
if file == 4
    f=ff;
end
fexcitation = f/1e6;
nof = length(fexcitation); % number of frequecies

% define gate for each freq
freq_len = total_time; %s time lapse for each freq
gate = 1:30000; % total trace length dont change
deltag = freq_len*fs;
gate_n = {};
% choose correct startpoint of the trace
startpoint =  4956;%5144%5026;
for iii = 1:nof
    gate_n{iii} = strcat(num2str(deltag*(iii-1)+startpoint),':',num2str(deltag*iii+startpoint-1));
end


%subfoldname = {'10kPa\','20kPa\','30kPa\','40kPa\','50kPa\','60kPa\','70kPa\','80kPa\',...
%'90kPa\','100kPa\',};


%% load and plot
Time=gate/fs;

color = {'r','b','k','m','y','g','r--','b--','k--','m--','y--','g'};
flag = 1;
for ii =1:length(subfoldname);
    % load water
    fp1 = [basepath,date,'/',foldname_1,'/',subfoldname{ii},'/'];
    dir_result = dir(fp1);
    load([fp1,dir_result(file).name]);
    wat = CompuScope_signals(2,:,gate);
    wat = squeeze(wat);
    
  
    % load bubble
    fp2 = [basepath,date,'/',foldname,'/',subfoldname{ii},'/'];
    dir_result = dir(fp2);
    load([fp2,dir_result(file).name]);
    bub = CompuScope_signals(2,:,gate);
    bub = squeeze(bub);
    
 if plotDetails == 1;
    figure(flag);
    plot(Time,wat); hold on;
    plot(Time,bub,'r');
    ylabel('voltage');
    xlabel('time [s]');
    legend('water','bubbles');
 end

    fs = 50e6;
    Nfft = length(CompuScope_signals);
    df = fs/Nfft;
    f = 0:df:fs-df;

    for nof_ii = 1:nof
        gate_p = str2num(gate_n{nof_ii});
        
  if plotDetails == 1;       
       figure(flag+1);
       subplot(3,7,nof_ii)
       plot(wat(gate_p));hold on
       plot(bub(gate_p),'r');
  end
                
            ftwat = abs(fft((wat(gate_p)),Nfft)).^2;
        
            ftbub = abs(fft((bub(gate_p)),Nfft)).^2;
                   
%             ftwat = abs(fft(detrend(wat(gate_p)),Nfft)).^2;
%         
%             ftbub = abs(fft(detrend(bub(gate_p)),Nfft)).^2;
%         
         

 if plotDetails == 1;
        %% plot attenuation coefficent
        figure(flag+2);
        plot(f./1e6,10*log10(ftwat/max(ftwat)),'b');hold on;
        plot(f./1e6,10*log10(ftbub/max(ftwat)),'LineWidth',4);
        xlim([fexcitation(1) fexcitation(end)]);hold on;
        xlabel('Frequency [MHz]');
        ylim([-20 0])
 end
        

        %% attenuation coefficent.
        fundamental_range = floor(fexcitation(nof_ii)*1e6/df-10):floor(fexcitation(nof_ii)*1e6/df+10);
        atten = ftwat./ftbub;
        acoeff(ii,nof_ii) = max(atten(fundamental_range));
    end
    atten_c(ii,:) = smooth(10*log10(acoeff(ii,:))/d);
   
     if plotDetails == 1;
    figure(flag+2)
    suptitle('Normalized Power Spectrum [dB]');
    flag = flag +4;
     end
end
%%%%%%%%%%%%%%%
figure
plotcolor(fexcitation,atten_c,subfoldname);
ylabel('Attenuation [dB/cm]')
xlabel('Frequency [MHz]')
xlim([0.7 3.3])
figformat

P = [2:2:8,10:10:80];
[attenA,maxAttenP] = max(atten_c');
frP = fexcitation(maxAttenP);
figure(20);
frP_fig = plot(P,frP);hold on;
xlabel('Pressure [kPa]')
ylabel('Fr(p) [MHz]')
title([date,foldname])
figformat

save(filename, 'atten_c', 'fexcitation');
