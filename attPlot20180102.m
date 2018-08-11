clear;%close all;
% ***************************************************
% cd attenuation fold. get attenuation data from water and bubble fold (cd each subfold of water and bubble fold)
% Plot of attenuation specturm for different excitation pressure (different attenuator db)
% Qian
%************************************************************

%% variables
basepath = '/Users/liqian/Google Drive/lab/lab/narrowband Attenuation/';
% fold and file info
date = '2016-08-05';
foldname_1 = 'wat'; % water
foldname = 'dppc5_1round'; % bubble
filename = [date,foldname];

gate = [2000:5000];
d =0.5; % cm
disp('change excitation info')
frange = [1:0.1:3];%MHz

% subfoldname1 = {'2kPa','3kPa','4kPa','5kPa','6kPa','7kPa','8kPa','9kPa','10kPa','20kPa','30kPa','40kPa','50kPa','60kPa','70kPa','80kPa'};%,'100kPa'}%,'30kPa\','30kPa\','30kPa\','30kPa\','30kPa\','30kPa\',};
%  subfoldname = {'2kPa','3kPa','4kPa','5kPa','6kPa','7kPa','8kPa','9kPa','10kPa','20kPa','30kPa','40kPa','50kPa','60kPa','70kPa','80kPa'};%,'100kPa'}%,'30kPa\','30kPa\','30kPa\','30kPa\','30kPa\','30kPa\',};
subfoldname1 = {'10kPa','20kPa','30kPa','40kPa','50kPa','60kPa','70kPa','80kPa'};%,'100kPa'}%,'30kPa\','30kPa\','30kPa\','30kPa\','30kPa\','30kPa\',};
 subfoldname = {'10kPa','20kPa','30kPa','40kPa','50kPa','60kPa','70kPa','80kPa'};%,'100kPa'}%,'30kPa\','30kPa\','30kPa\','30kPa\','30kPa\','30kPa\',};

 color = {'r-','b-','k-','m-','g-','r--','b--','k--','m--','g--','r-*','b-*','k-*','g-*','m-*','y-*','r--','b--','k--','m--','y--','g'};


%subfoldname1 ={'10kPa','30kPa','50kPa','70kPa','90kPa'};
 %subfoldname ={'10kPa','30kPa','50kPa','70kPa','90kPa'};
% define gate for each freq
% gate = 4000:30000; % total trace length

%subfoldname = {'10kPa\','20kPa\','30kPa\','40kPa\','50kPa\','60kPa\','70kPa\','80kPa\',...
%'90kPa\','100kPa\',};

%% load and plot

flag = 1;
for ii =1:length(subfoldname) % pressure fold
    % load water at different f same p
    fp1 = [basepath,date,'/',foldname_1,'/',subfoldname1{ii},'/'];
    dir_result = dir(fp1);
    
    for kk = 3:length(dir_result)
    load([fp1,dir_result(kk).name]);
    wat = Y;
    

%     figure(flag);
%     plot(wat); hold on;
    % %%%%%%%%%%%%%%
    % sig_w  = reshape(wat(1,:),gate_len,nof);
    % %for ii = 1:nof;figure;plot(sig_w(:,1));end;
    % %%%%%%%%%%%%%%%%
    % load bubble
    fp2 = [basepath,date,'/',foldname,'/',subfoldname{ii},'/'];
    dir_result1 = dir(fp2);
    load([fp2,dir_result1(kk).name]);
     bub = Y;
    

%     figure(flag);
%     plot(bub,'r');
%     ylabel('voltage');
%     xlabel('time [s]');
%     legend('water','bubbles');
fs = 1/(T(2)-T(1));
Nfft = 1e4;
df = fs/Nfft;
f = 0:df:fs-df;


            ftwat = abs(fft(detrend(wat),Nfft)).^2;      
            %figure(1);plot(f,ftwat,color{kk-2});hold on;
        
            ftbub = abs(fft(detrend(bub),Nfft)).^2;
            %figure(1);plot(f,ftbub,color{kk-2});hold on;
        

%                 plot(wat(1,gate_p));hold on
%                 plot(bub(1,gate_p),'r');


        %% plot attenuation coefficent

%      figure(100+ii);
        fundamental_range = floor((frange(kk-2)-0.1)*1e6/df):floor((frange(kk-2)+0.1)*1e6/df);
        watpower = 10*log10(ftwat/max(ftwat(fundamental_range)));
        bubpower = 10*log10(ftbub/max(ftwat(fundamental_range)));
%         plot(f(fundamental_range)./1e6,watpower(fundamental_range),'b');hold on;
%         plot(f(fundamental_range)./1e6,bubpower(fundamental_range),'r');
%         xlim([0 5]);
%         hold on;
%         xlabel('Frequency [MHz]');
%         ylim([-20 0])
        
        %legend('water','bubbles');
        %         figure(flag+3);
        %         subplot(5,2,nof_ii)
        %         plot(f./1e6,10*log10(ftwat_mean./ftbub_mean)/d,color{ii});xlim([0 3]);hold on
        %         xlabel('Frequency [MHz]');


        %% attenuation coefficent.
        acoeff(ii,kk-2) = max(watpower(fundamental_range))-max(bubpower(fundamental_range));
        
        flag = flag +3;
        
    end
    %         figure(flag+3)
    %         suptitle('Attenuation [dB/m]');
    atten_c(ii,:) = smooth(acoeff(ii,:))/d;
    fexcitation = frange;
end
%%%%%%%%%%%%%%%
plotcolor(fexcitation,atten_c,subfoldname);
ylabel('Attenuation [dB/cm]')
xlabel('Frequency [MHz]')
figformat

save(filename, 'atten_c', 'fexcitation');