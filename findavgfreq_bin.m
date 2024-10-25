function [freqvecmean, allfreq] = findavgfreq_bin(timecourse, threshforplot, binsize)

%timecourse = timecourse(500:end,:);
s = timecourse; %use calcium to determine phase
t = (0:size(s,1)-1)*binsize/10; %time vector
Ts = mean(diff(t)); %% Sampling period
Fs = 1/Ts;  % Sampling Frequency
Fn = Fs/2;  %Nyquist Frequency
L  = size(s,1); %length of signal
cellnumber = size(s,2);

fts = fft(s-(ones(L,1)*mean(s,1)));
%fts1 = fft(s-mean(s,1));

for i=1:size(timecourse,2)
    ft_ind = fts(:,i);
    course =s(:,i)-mean(s(:,i));
    f = Fs*(0:(ceil(L/2)))/L;
    
    ft_test = fft(s(:,i)-mean(s(:,i)));
    if (ft_test~=ft_ind)
    end
    
    P2 = abs(ft_ind/L);
    P1 = P2(1:(ceil(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);
    
    %     figure
    %     plot(f, P1)
    %     xlim([0 1])
    
    amp_fts = abs(ft_ind);      % Spectrum Amplitude
    phs_fts = angle(ft_ind);
    highestfreq = f(find(P1==max(P1)));%this is the prob
    value = find(f==highestfreq);
    freqvec(i) = f(value);
    ampvec(i) = amp_fts(value);
    phasevec(i) = unwrap(phs_fts(value));
    
    sortedpower = sort(P1,'descend');
    
    %%plot the frequencies back on time course as sin wave
    %         figure
    %         plot(s(:,i)-mean(s(:,i)))
    %          hold on
    %         plot(max(P1)*sin(2*pi*highestfreq *t))
    %         hold on
    %         hold on
    %         f_2nd = f(find(P1== sortedpower(2)));
    %         plot(sortedpower(2)*sin(2*pi*f_2nd *t))
    %         hold on
    %         f_3rd = f(find(P1== sortedpower(3)));
    %         plot(sortedpower(3)*sin(2*pi*f_3rd *t))
    
    %     figure()
    %     plot(timecourse(:,i))
    %     [pks,locs,w,p]  = findpeaks(course);
    %     peaks = find(w>50);
    %     if length(peaks)<2
    %         if length(p) > 2
    %             peaks = find(w<50);
    %             if length(peaks) >2
    %                 peaks = [peaks(length(peaks)/2), peaks(length(peaks))];
    %             end
    %             manualfreq(i) = 1/(t(locs(peaks(2))) - t(locs(peaks(1))));
    %
    %             figure
    %             plot(s(:,i)-mean(s(:,i)))
    %             hold on
    %             plot(locs(peaks),pks(peaks), 'ro');
    %             hold on
    %             plot(max(P1)*sin(2*pi*manualfreq(i) *t))
    %             hold on
    %         else
    
    %             figure
    %             plot(s(:,i)-mean(s(:,i)))
    %             hold on
    %     try
    plotvalues = [1 :50:1000];
    if max(s(:,i)-mean(s(:,i)))> .1*10^-4

            [pks,locs,w,p]  = findpeaks(course);
            maxcourse = max(course);
            %findtops = abs(maxcourse-pks);
            peaks = 1:length(locs);
            %peaks = find(w(peaks)>10);
            
            highvalues = zeros(1, length(course));
            midwaypoint = zeros(1, length(course));
            halfwaypoint = min(course) + (max(course)-min(course))/2;
            midwaypoint(find(course>halfwaypoint))=1;
            middiff = abs(diff(midwaypoint));
            wheremid = find(middiff==1);
            minmid = wheremid(1);
            minni=2;
            while minmid < locs(peaks(1))
                minmid = wheremid(minni);
                minni=minni+1;
            end
            
            thecorrectvalue = find(locs(peaks)>minmid);
            
            if length(peaks)>2
                thecorrectvalue = find(locs(peaks)>minmid);
                newpeaks = [peaks(1),peaks(thecorrectvalue(1))];
                more=2;
                while (locs(newpeaks(2))-locs(newpeaks(1))) < 20/(binsize)
                    newpeaks = [peaks(1),peaks(thecorrectvalue(more))];
                    more=more+1;
                    if more > length(thecorrectvalue)
                        'break';
                        break
                    end
                end
                peaks=newpeaks;
                newfreq = 1/(t(locs(peaks(2))) - t(locs(peaks(1))));
                    manualfreq(i)=newfreq;
                %                 figure
                %                 plot(s(:,i)-mean(s(:,i)))
                %                 hold on
                %                 plot(locs(peaks),pks(peaks), 'ro');
                %                 hold on
                %                 plot(max(P1)*sin(2*pi*manualfreq(i) *t))
                %                 title(['new manualfreq = ' num2str(newfreq)]);
            end
                if ismember(i, plotvalues)
                    figure
                    plot(s(:,i)-mean(s(:,i)));
                    hold on
                    plot(locs(1:2),pks(1:2), 'ro');
                    hold on
                    plot(max(P1)*sin(2*pi*manualfreq(i) *t));
                    hold on
                    disp(['high freq i=' num2str(i)]);
                    title(['orig manualfreq = ' num2str(manualfreq(i))]);
                end
    else
        manualfreq(i) = nan;
        %             figure
        %             plot(s(:,i)-mean(s(:,i)))
        %             disp(['freq 0 i=' num2str(i)]);
    end
    
end

freqvecmean= nanmean(manualfreq);
allfreq = manualfreq;

end

