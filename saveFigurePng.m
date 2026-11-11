function saveFigurePng(fig, pathNoExt, dpi)
%SAVEFIGUREPNG  Save a figure as a high-resolution PNG (best format for PowerPoint: lossless, small, sharp text).
%   saveFigurePng(fig, '/path/name')            -> /path/name.png  at 300 dpi
%   saveFigurePng(fig, '/path/name', 600)       -> 600 dpi
if nargin < 3, dpi = 300; end
[p,n,~] = fileparts(pathNoExt);
f = fullfile(p,[n '.png']);
try
    exportgraphics(fig,f,'Resolution',dpi,'BackgroundColor','white');
catch
    try
        print(fig,f,'-dpng',sprintf('-r%d',dpi));
    catch ME
        warning('Could not save %s (%s).',f,ME.message);
    end
end
end
