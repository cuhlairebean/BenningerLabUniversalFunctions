function dutycycle = getDC(datashort,threshholdDC)

totaltime = size(datashort,1);

figure
plot(datashort)
hold on
plot(ones(totaltime,1)*threshholdDC)
title('Data vs Duty Cycle Threshhold')

DCvec = zeros(size(datashort));
DCvec(find(datashort > threshholdDC)) = 1;
amt_on = sum(DCvec); 

dutycycle = mean(amt_on)/totaltime;

end

