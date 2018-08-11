%% use for mb data <= 3mhz. not only for c20
%% work with test.m. plot exp sh data vs pressure amd freq. minus noise
    for k = 1:11%length(excitation_f) %excitation freq. one freq,multiple pressure curve
        figure(k);
        lh=findall(gca,'type','line');
        xx=get(lh,'xdata');
        yy=get(lh,'ydata');
        yyy = zeros(length(yy)/2,500);
        
        for Pcurveii = 1:length(yy)/2 % pressure curve
            % minus 6um standard beads
            new = yy{Pcurveii,:}-yy{Pcurveii+length(yy)/2,:};
            yyy(Pcurveii,1:length(new)) = new;
            SHf_gate = floor([(0.5*excitation_f(k)-0.2)*1e6:df:(0.5*excitation_f(k)+0.2)*1e6]/df+2);
            subharmonicdB(k,Pcurveii) = sqrt(10^(max(yyy(Pcurveii,SHf_gate))/10));
            figure(20+k)
            plot(xx{Pcurveii,1}/frequency(k),yyy(Pcurveii,1:length(new)));
            hold on;
        end
       xlim([0.2 1.2]);
        %ylim([-30 0])
        figformat;
        grid on
        xlabel('Frequency/Fexcitation')
        ylabel('Power spectrum (dB)')
        legend(strcat(cellfun(@num2str,num2cell(PP(end:-1:1)), 'UniformOutput', false),'kPa'))
        title(['f = ',sprintf('%.2f',frequency(k)),'minus noise'])
    end
 
%% subharmonic amplitude |A|; ignor dB in the subharmonicdB name
figure(107);% 
for ii = 1:length(subharmonicdB) % freq
        plot(20*log10(fliplr(subharmonicdB(ii,:))));
        %colormap(jet)
        %xlabel('Frequency (MHz)')
        xlabel('Excitation Pressure (kPa)')
        title(['subharmoinc response - not normalized with fundamental response',mark]);
        hold on;
end
figformat
legend('1MHz','1.2MHz','1.4MHz','1.6MHz','1.8MHz','2MHz','2.2MHz','2.4MHz','2.6MHz','2.8MHz','3MHz')

 
 %% normalize with maximum f p RESPONSE
 
        figure(107);
        lh=findall(gca,'type','line');
        xx=get(lh,'xdata');
        yy=get(lh,'ydata');
        
        for Pcurveii = 1:length(yy) % pressure curve 
            ny(Pcurveii,:) = yy{Pcurveii,:}/max([yy{:}]);
            figure(108)
            plot(xx{Pcurveii,1}.*10,ny(Pcurveii,:));hold on
        end
        figformat
ylabel('sub |A| normalized with maximum response')
xlabel('pressure [kpa]')
legend('3Mhz','2.8Mhz','2.6Mhz','2.4Mhz','2.2Mhz','2.0Mhz','1.8Mhz','1.6Mhz','1.4Mhz','1.2Mhz','1Mhz')        
 