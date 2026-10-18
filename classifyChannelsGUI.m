function [gCh, rCh, labels] = classifyChannelsGUI(stack, fileLabel)
% CLASSIFYCHANNELSGUI  Show every channel of an HxWxC stack and let the user label each one.
%   [gCh, rCh, labels] = classifyChannelsGUI(stack, fileLabel)
%   Each channel is displayed (all side by side, in its own colour-tinted view) and you are
%   asked, channel by channel, whether it is Green, Red or Other. Exactly one Green and one Red
%   channel must be chosen (you are asked again otherwise).
%   gCh/rCh : channel numbers (NaN if cancelled);  labels : 1xC cell of 'Green'/'Red'/'Other'
if nargin < 2, fileLabel = ''; end
nC = size(stack,3);
gCh = NaN; rCh = NaN; labels = repmat({'Other'},1,nC);
opts = {'Green','Red','Other (ignore)'};

fig = figure('Name',['Channels: ' fileLabel],'NumberTitle','off', ...
    'Units','normalized','Position',[0.05 0.2 0.9 0.55]);
ax = gobjects(1,nC);
for c = 1:nC
    ax(c) = subplot(1,nC,c,'Parent',fig);
    im = double(stack(:,:,c));
    im = (im - min(im(:)))/max(eps,(max(im(:))-min(im(:))));
    imshow(im,'Parent',ax(c));
    title(ax(c),sprintf('Channel %d',c));
end

while true
    for c = 1:nC
        set(ax(c),'Visible','on','XColor','y','YColor','y','LineWidth',3);
        for o = setdiff(1:nC,c), set(ax(o),'XColor','none','YColor','none'); end
        title(ax(c),sprintf('Channel %d  <-- classify',c),'Color','r');
        drawnow;
        sel = listdlg('PromptString',sprintf('Which colour is Channel %d?',c), ...
            'ListString',opts,'SelectionMode','single','InitialValue',min(c,3), ...
            'Name','Classify channel','ListSize',[220 80]);
        if isempty(sel)
            if isvalid(fig), close(fig); end
            labels = repmat({'Other'},1,nC); return;   % cancelled
        end
        title(ax(c),sprintf('Channel %d: %s',c,strtok(opts{sel})),'Color','k');
        labels{c} = strtok(opts{sel});
    end
    nG = nnz(strcmp(labels,'Green')); nR = nnz(strcmp(labels,'Red'));
    if nG == 1 && nR == 1
        gCh = find(strcmp(labels,'Green')); rCh = find(strcmp(labels,'Red'));
        break;
    end
    uiwait(msgbox(sprintf('Exactly ONE Green and ONE Red channel are required (you chose %d Green, %d Red). Please classify again.',nG,nR),'Classify channels','modal'));
end
if isvalid(fig), close(fig); end
end
