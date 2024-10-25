function [maxcorrarea, Dat] = getcoordarea_forJosh(ActMask, LogicalMap, images,T)
[sx, sy] = size(ActMask);
CorrPeakAmp = zeros(sx,sy);
TopCorrPeakAmp = zeros(sx,sy);
    if (max(max(ActMask))>0)
        [Dat]=HIanalyzeU(images,T,0,ActMask);
        Dat.images=images;
        Dat.ActMask = ActMask;
        Dat.T=T;
        Dat=CoordinatedAreas_forJosh(Dat);
        
        CorrPeakAmp = Dat.Img();
        %CorrPeakAmp = CorrPeakAmp + ActMask;
        CorrPeakAmp(~LogicalMap)=0;
        
        STATS=regionprops(CorrPeakAmp,'area');
        for j=1:length(STATS)
            stats(j)=STATS(j).Area;
        end
        statsordered = sort(stats, 'descend');
        
        maxstatsarea = statsordered(1);
        len = length(statsordered);
        if len < 5
            value = len;
        else
            value = 5;
        end
%             r = find(stats>=statsordered(value));
%             gettopstat = find(stats>=statsordered(1));
%             number = 1;
%             for i=r
%                 TopCorrPeakAmp(CorrPeakAmp==i)=number;
%                 if i==gettopstat
%                     Dat.maxcorrarea_num = number;
%                 end
%                 number = number+1;
%             end
            
            %% order them by size
            count = 1;
            for i=1:length(statsordered)
                try
                r = find(stats==statsordered(i));
                TopCorrPeakAmp(CorrPeakAmp==r)=count;
                count = count+1;
                catch
                end
                
            end
            
        maxcorrarea = max(stats);
    else
        maxcorrarea = 0;
        Dat.AMCA = 0;
    end
    figure()
    imagesc(CorrPeakAmp);
    colorbar;
    Dat.CorrPeakAmp = CorrPeakAmp;
%     
%     
    figure()
    imagesc(TopCorrPeakAmp);
    colorbar;
    Dat.TopCorrPeakAmp = TopCorrPeakAmp;
    

end

