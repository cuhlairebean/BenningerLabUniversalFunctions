function [ dataOut ] = GetPhasev2(imagesfull, masks,zstacks, PixelSize, second, starttime, endtime);
%HIANALYZE Summary of this function goes here
%   Detailed explanation goes here
%images=cat(3,images{:});
if starttime == -1
    st=1;
else
    st = starttime;
end

if endtime == -1
    ed=size(imagesfull,3);
else
    ed=endtime;
end

for z=1:zstacks
    
    images=double(imagesfull(:,:,:,z));
    
    IsletAvg=mean(mean(images,1),2);
    IsletAvg = IsletAvg(:);
    
    figure; plot(IsletAvg/mean(IsletAvg));
    
    images=images(:,:,st:ed);
    mask = masks(:,:,z);
    imagesRaw=mat2gray(images);
    % matrixClusterLo=imagesRaw;
    % h = fspecial('average', 1);
    % avgIslet=imfilter(matrixClusterLo,h,'replicate');
    % images=avgIslet;
    
    imageMean=mean(images,3);
    img=images(:,:,1);
    sx=size(img,1);
    sy=size(img,2);
    sz=size(images,3);
    
    corrMap = nan(sx,sy);
    phaseMap = nan(sx, sy);
    hilbertMap =  nan(sx, sy);
    
    for i=1:sz
        images(:,:,i)=images(:,:,i).*mask;
        %images(:,:,i)=imfilter(images(:,:,i),[5 5]).*mask;%%%%%%%%%%%%%%%%%%%%%%
    end
    
    DM=mean(images,1);
    DM=mean(DM,2);
    DMean=mean(mean(images,1),2);
    DMean = DMean(:);
    
    figure; plot(DMean/mean(DMean));
    hold on
 
    timeVec = 1:second:length(DMean)*second; %time vector
    Ts = mean(diff(timeVec)); %% Sampling period
    Fs = 1/Ts;  % Sampling Frequency
    L = sz;
    NFFT=2^nextpow2(L);
    TF=Fs./2*linspace(0,1,NFFT/2+1);
    TF(1)=[];
    TF=1./TF;
    Fn = Fs/2;  %Nyquist Frequency
    
    DF=fft(DMean-mean(DMean));
   
    f = Fs*(0:(ceil(L/2)))/L;
    P2 = abs(DF/L);
    P1 = P2(1:(ceil(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);
    
    P1ordered = sort(P1);
    highestfreq = f(find(P1==max(P1)));%this is the prob
%     value = find(f==highestfreq);
%     amp_fts = abs(DF);      % Spectrum Amplitude
%     phs_fts = angle(DF);
%     freqvec(i) = f(value);
%     ampvec(i) = amp_fts(value);
%     phasevec(i) = unwrap(phs_fts(value));
   
    [x,y]=find(mask==1);
    
    for i=1:length(x)
        x1 = x(i)-10;
         x2 = x(i)+10;
          y1 = y(i)-10;
           y2 = y(i)+10;
        
        if x1 < 1
            x1=1;
        end
         if x2 > sx
            x2=sx;
         end
         if y1 <1
            y1=1;
         end
         if y2 > sy
            y2=sy;
        end
        
        dat=images(x1:x2,y1:y2,:);
        dat=mean(dat,2);
        dat=mean(dat,1);
        dat=dat(:);
        
        %dat=detrend(dat);
        Y=fft(dat-mean(dat));
        
        P2 = abs(Y/L);
        P1 = P2(1:(ceil(L/2)+1));
        P1(2:end-1) = 2*P1(2:end-1);
         
        amp_fts = abs(Y);      % Spectrum Amplitude
        phs_fts = angle(Y);
        value1 = find(f==highestfreq);
%         freqvec(i) = f(value1);
%         ampvec(i) = amp_fts(value1);
%         phasevec(i) = unwrap(phs_fts(value1));
%         phasevec1(i) = angle(Y(value1));
    
        R = corrcoef(dat,DMean);
        if R(2) >.8
            if rand(1) >.9
                plot(dat/mean(dat));
                hold on
            end
            
            [c lags]=xcov(dat,DMean,'coeff');
            [maxCV maxCL]=max(c);
            
            if ~isnan(dat(1))
                corrMap(x(i),y(i))=lags(maxCL);
            end
            phaseMap(x(i),y(i))=angle(Y(value1));
            %hilbertMap(x(i),y(i)) = unwrap(angle(hilbert(detrend(dat))));
        end
    end
    phaseMap=phaseMap/(2*pi*highestfreq);
    
    maxp = max(max(phaseMap));
    minp = min(min(phaseMap));
    tenp = (maxp-minp)*.1;
    
    [row col] = find(phaseMap > maxp-tenp);
    [row2 col2] = find(phaseMap < minp+tenp);
    
    count = 1;
    for ii = 1:length(row2)
        for jj = 1:length(row)
            distance(count) = norm([row(jj) col(jj)] - [row2(ii) col2(ii)]);
            maxpall(count) = phaseMap(row(jj), col(jj));
            minpall(count) = phaseMap(row2(ii), col2(ii));
            count = count+1;
        end
    end
    distance = distance.*PixelSize;
    
    allspeed = distance./(maxpall-minpall);
    dataOut.timediff(z) = maxp-minp;
    dataOut.frequency(z) = highestfreq;
    
    dataOut.meanspeed(z) = mean(allspeed);
    dataOut.medspeed(z)  = median(allspeed);
    dataOut.maxspeed(z) = max(allspeed);
    dataOut.minspeed(z) = min(allspeed);
    
    dataOut.meandis(z) = mean(distance);
    dataOut.meddis(z)  = median(distance);
    dataOut.maxdis(z) = max(distance);
    dataOut.mindis(z) = min(distance)
    
    %dataOut.corrMap(:,:,z)=corrMap;
    
    figure; 
    imagesc(corrMap);
    colorbar
    title('corrMap')
    
    %dataOut.phaseMap(:,:,z)=phaseMap;
    
    figure; 
    imagesc(phaseMap);
    colorbar
    title('phaseMap')
    
end

end

