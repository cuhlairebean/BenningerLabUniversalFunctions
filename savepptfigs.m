function [] = savepptfigs(fullpath, filename, fulltimeboolean)

if fulltimeboolean
    todaysdate = datestr(datetime('now'));
else
    todaysdate = datestr(datetime('today'));
end
todaysdate = regexprep(todaysdate,'-','');
todaysdate = regexprep(todaysdate,' ','_');
todaysdate = regexprep(todaysdate,':','');
todaysdate

pptname = [fullpath '\' filename '_' todaysdate 'ppt'];
saveAllFigsToPPT(pptname)
disp(['Figures saved ' fullpathname]);

%save([Path 'workspace' datestr(datetime('today')) '.mat'])
%save(savedfile, 'variable');
%save([savedfile '_Data.mat'], 'filename');
%save([Path 'workspace' datestr(datetime('today')) '.mat'])


end

