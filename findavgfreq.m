function [freqvecmean, allfreq] = findavgfreq(timecourse, threshforplot)

timecourse = timecourse(1:end,:);
s = timecourse; %use calcium to determine phase
t = (0:size(s,1)-1)/100; %time vector
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
    
   
    plotvalues = [1 7 10 230 687 948];
    if max(s(:,i)-mean(s(:,i)))> .1*10^-4
        highvalues = zeros(1, length(course));
        midwaypoint = zeros(1, length(course));
        halfwaypoint = min(course) + (max(course)-min(course))/2;
        midwaypoint(find(course>halfwaypoint))=1;
        middiff = abs(diff(midwaypoint));
        wheremid = find(middiff==1);
        minmid = wheremid(1);
        minni=2;
        while minmid < 20
            minmid = wheremid(minni);
            minni=minni+1;
        end
        if course(2)- course(1) < 0
            highvalues(find(course < course(1)))=1;
            findchange = find(diff(highvalues)==1);
        else
            highvalues(find(course > course(1)))=1;
            findchange = find(diff(highvalues)==1);
        end
        
        if length(findchange)>2
            thecorrectvalue = find(findchange>minmid);
            try
            locs = [findchange(1),findchange(thecorrectvalue(1))];
            more=2;
            while (locs(2)-locs(1)) < 20
                locs = [findchange(1),findchange(thecorrectvalue(more))];
                more=more+1;
                if more > length(thecorrectvalue)
                    'break';
                    break
                end
            end
            catch 
                figure
            plot(s(:,i)-mean(s(:,i)));
            disp(['error on i=' num2str(i)]);
            end
        elseif length(findchange)==1
            highvalues = zeros(1, length(course));
            highvalues(find(course <halfwaypoint))=1;
            findchange = find(diff(highvalues)==1);
            locs = findchange;
            if length(findchange)==1
                highvalues = zeros(1, length(course));
                highvalues(find(course <halfwaypoint))=1;
                findchange = find(diff(highvalues)==-1);
                locs = findchange;
            elseif length(findchange)>2
                thecorrectvalue = find(findchange>minmid);
                locs = [findchange(1),findchange(thecorrectvalue(1))];
                more=2;
                while (locs(2)-locs(1)) < 20
                    locs = [findchange(1),findchange(thecorrectvalue(more))];
                    more=more+1;
                end
            end
            %                 plot(s(:,i)-mean(s(:,i)))
            %                 hold on
            %                 plot(findchange,[halfwaypoint, halfwaypoint], 'ro');
            %                 hold on
            %                 plot(max(P1)*sin(2*pi*manualfreq(i) *t))
        else
            locs = findchange;
        end
        try
            manualfreq(i) = 1/(t(locs(2)) - t(locs(1)));
            pks = [course(1), course(1)];
        catch
             manualfreq(i) = nan;
                    figure
                    plot(s(:,i)-mean(s(:,i)))
                    disp(['freq nan i i=' num2str(i)]);
            
        end
        
        %plot funky ones
        %                        plot(locs,pks, 'ro');
        %                        hold on
        %                        plot(max(P1)*sin(2*pi*manualfreq(i) *t))
        %                        hold on
        
        if manualfreq(i) < .005
            figure
            plot(s(:,i)-mean(s(:,i)));
            disp(['low freq i=' num2str(i)]);
        end
        %outputgraphs = 1:40:1000;
        if manualfreq(i) > threshforplot %%for continuous
            %if manualfreq(i) > .21 %%for bimodal
%             figure
%             plot(s(:,i)-mean(s(:,i)));
%             hold on
%             plot(locs(1:2),pks, 'ro');
%             hold on
%             plot(max(P1)*sin(2*pi*manualfreq(i) *t));
%             hold on
%             disp(['high freq i=' num2str(i)]);
%             title(['orig manualfreq = ' num2str(manualfreq(i))]);
            
            %plot(timecourse(:,i))
            %try
            [pks,locs,w,p]  = findpeaks(course);
            maxcourse = max(course);
            findtops = abs(maxcourse-pks);
            peaks = find(findtops<.3*10^-4);
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
                while (locs(newpeaks(2))-locs(newpeaks(1))) < 20
                    newpeaks = [peaks(1),peaks(thecorrectvalue(more))];
                    more=more+1;
                    if more > length(thecorrectvalue)
                        'break';
                        break
                    end
                end
                peaks=newpeaks;
                newfreq = 1/(t(locs(peaks(2))) - t(locs(peaks(1))));
                if newfreq < manualfreq(i)
                    manualfreq(i)=newfreq;
                end
%                 figure
%                 plot(s(:,i)-mean(s(:,i)))
%                 hold on
%                 plot(locs(peaks),pks(peaks), 'ro');
%                 hold on
%                 plot(max(P1)*sin(2*pi*manualfreq(i) *t))
%                 title(['new manualfreq = ' num2str(newfreq)]);
            end
        end
%         if ismember(i, plotvalues)
%             figure
%             plot(s(:,i)-mean(s(:,i)));
%             hold on
%             plot(locs(1:2),pks(1:2), 'ro');
%             hold on
%             plot(max(P1)*sin(2*pi*manualfreq(i) *t));
%             hold on
%             disp(['high freq i=' num2str(i)]);
%             title(['orig manualfreq = ' num2str(manualfreq(i))]);
%         end
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

