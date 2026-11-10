function [roiLabel, info] = drawROIsGUI(baseImg, varargin)
%DRAWROISGUI  Simple keyboard-driven GUI to draw many freehand ROIs on an image.
%
%   [roiLabel, info] = drawROIsGUI(baseImg, 'Name','Cell', 'Overlay',stack, 'ParentMask',IsletMask)
%
%   No on-screen buttons: everything is a key press (click the figure once so it
%   has keyboard focus; press 'h' for the key list).
%
%   KEYS
%     =  /  -        zoom in / out          arrow keys   pan          r  reset view
%     b  /  B        base image brighter / dimmer
%     n  /  N        base image more / less contrast
%     o  /  O        overlay brighter / dimmer
%     p  /  P        overlay more / less contrast
%     t  /  T        overlay more / less opaque (translucent)
%     1 2 3 4 ...    toggle overlay channel 1,2,3,4...      v  overlay on/off
%     d              draw next ROI (click-drag freehand, release to finish)
%     u              undo last ROI                          Enter / f   FINISH
%
%   The overlay (any resolution, e.g. a high-res _pre.czi) is stretched onto the
%   base image's coordinate system, so zooming/panning can never misalign it.
%
%   INPUTS (name/value): 'Name','Title','Colormap','Overlay' (HxWxC),
%       'ChannelNames','ParentMask' (label matrix),'ParentName'
%   OUTPUTS
%     roiLabel  label matrix, size(baseImg): pixel = ROI number in draw order
%               (later ROIs overwrite earlier ones where they overlap)
%     info      .masks (cell of logical masks), .parentIdx (Nx1; parent # each ROI
%               falls in, 0 = none), .parentOverlap (Nx1 fraction inside parent)

p = inputParser;
addParameter(p,'Name','ROI');
addParameter(p,'Title','');
addParameter(p,'Colormap','gray');
addParameter(p,'Overlay',[]);
addParameter(p,'ChannelNames',{});
addParameter(p,'ParentMask',[]);
addParameter(p,'ParentName','Islet');
parse(p,varargin{:});
opt = p.Results;

baseImg = double(baseImg);
[sx, sy] = size(baseImg);
bmin = min(baseImg(:)); bmax = max(baseImg(:));
if bmax <= bmin, bmax = bmin + 1; end
bCenter = (bmin+bmax)/2; bWidth = bmax-bmin;      % base display window

chColors = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 0];
chNames  = {'red','green','blue','magenta','cyan','yellow'};

masks = {}; hLines = {}; hTexts = {}; parentIdx = zeros(0,1); parentFrac = zeros(0,1);
chOn = []; ovNorm = {}; oCenter = 0.5; oWidth = 1; oAlpha = 0.5; ovVisible = true;
status = '';

hasParent = ~isempty(opt.ParentMask);
nParents = 0;
if hasParent
    opt.ParentMask = double(opt.ParentMask);
    nParents = max(opt.ParentMask(:));
end
pcols = hsv(max(nParents,1));

%% figure
fig = figure('Name',[opt.Name ' ROI GUI'],'NumberTitle','off', ...
    'Units','normalized','Position',[0.1 0.08 0.7 0.85]);
ax = axes('Parent',fig,'Position',[0.04 0.05 0.68 0.85]);
keyMenu = { ...
    'ACTION KEYS', ...
    '', ...
    'VIEW', ...
    '  =  /  -     zoom in / out', ...
    '  arrows      pan', ...
    '  r           reset view', ...
    '', ...
    'BASE IMAGE', ...
    '  b / B       brighter / dimmer', ...
    '  n / N       more / less contrast', ...
    '', ...
    'OVERLAY (_pre.czi)', ...
    '  o / O       brighter / dimmer', ...
    '  p / P       more / less contrast', ...
    '  t / T       more / less opaque', ...
    '  1 2 3 4     toggle channel', ...
    '  v           overlay on / off', ...
    '', ...
    'ROIs', ...
    '  d           draw next ROI', ...
    '  u           undo last ROI', ...
    '  Enter / f   FINISH', ...
    '', ...
    '  h           print keys in', ...
    '              Command Window'};
uicontrol(fig,'Style','text','Units','normalized','Position',[0.74 0.05 0.25 0.85], ...
    'String',keyMenu,'HorizontalAlignment','left','FontName','Courier', ...
    'FontSize',11,'BackgroundColor',get(fig,'Color'));

hBase = imagesc(ax,[1 sy],[1 sx],baseImg);
colormap(ax,opt.Colormap);
axis(ax,'image'); hold(ax,'on');
set(ax,'XLim',[0.5 sy+0.5],'YLim',[0.5 sx+0.5]);
hOv = image(ax,'XData',[1 sy],'YData',[1 sx],'CData',zeros(2,2,3),'AlphaData',0,'Visible','off');

if ~isempty(opt.Overlay), setOverlay(opt.Overlay); end

% parent outlines
if hasParent
    for k = 1:nParents
        [r,c] = find(opt.ParentMask==k);
        if isempty(r), continue; end
        B = bwboundaries(opt.ParentMask==k);
        for b = 1:numel(B)
            plot(ax,B{b}(:,2),B{b}(:,1),'--','Color',pcols(k,:),'LineWidth',1.2);
        end
        text(ax,mean(c),mean(r),sprintf('%s %d',opt.ParentName,k),'Color',pcols(k,:), ...
            'FontWeight','bold','HorizontalAlignment','center');
    end
end

setappdata(fig,'done',false);
set(fig,'WindowKeyPressFcn',@onKey,'CloseRequestFcn',@onClose);
updateBase(); showHelp(); setStatus('Press h for keys. d = draw ROI, Enter = finish.');

% wait until finished
while ~getappdata(fig,'done')
    pause(0.05);
    if ~isvalid(fig), break; end
end

roiLabel = zeros(sx,sy);
for n = 1:numel(masks), roiLabel(masks{n}) = n; end
info.masks = masks; info.parentIdx = parentIdx; info.parentOverlap = parentFrac;
if isvalid(fig), delete(fig); end

%% ---------------- local (nested) functions ----------------
    function onKey(~,e)
        ch = e.Character; k = e.Key;
        if isempty(ch) && isempty(k), return; end
        switch k
            case 'leftarrow',  panView(-0.2,0);
            case 'rightarrow', panView(0.2,0);
            case 'uparrow',    panView(0,-0.2);
            case 'downarrow',  panView(0,0.2);
            case {'return','enter'}, finish();
            otherwise
                switch ch
                    case {'=','+'}, zoomView(1/1.5);
                    case {'-','_'}, zoomView(1.5);
                    case 'r', set(ax,'XLim',[0.5 sy+0.5],'YLim',[0.5 sx+0.5]);
                    case 'b', bCenter = bCenter - 0.1*bWidth; updateBase();
                    case 'B', bCenter = bCenter + 0.1*bWidth; updateBase();
                    case 'n', bWidth = bWidth/1.25; updateBase();
                    case 'N', bWidth = bWidth*1.25; updateBase();
                    case 'o', oCenter = oCenter - 0.1*oWidth; updateOverlay();
                    case 'O', oCenter = oCenter + 0.1*oWidth; updateOverlay();
                    case 'p', oWidth = oWidth/1.25; updateOverlay();
                    case 'P', oWidth = oWidth*1.25; updateOverlay();
                    case 't', oAlpha = min(1,oAlpha+0.1); updateOverlay();
                    case 'T', oAlpha = max(0,oAlpha-0.1); updateOverlay();
                    case 'v', ovVisible = ~ovVisible; updateOverlay();
                    case 'u', undoROI();
                    case 'd', addROI();
                    case 'f', finish();
                    case 'h', showHelp();
                    otherwise
                        n = str2double(ch);
                        if ~isnan(n) && n >= 1 && n <= numel(chOn)
                            chOn(n) = ~chOn(n); updateOverlay();
                        end
                end
        end
    end

    function updateBase()
        lo = bCenter - bWidth/2; hi = bCenter + bWidth/2;
        if hi <= lo, hi = lo + 1e-6; end
        set(ax,'CLim',[lo hi]);
    end

    function setOverlay(stack)
        stack = single(stack);
        [H,W,nC] = size(stack);
        ovNorm = cell(1,nC);
        for c = 1:nC
            chn = stack(:,:,c);
            v = sort(chn(1:max(1,floor(numel(chn)/100000)):end));
            top = v(max(1,round(0.998*numel(v)))); lo = v(1);
            if top <= lo, top = lo + 1; end
            ovNorm{c} = min(max((chn-lo)./(top-lo),0),1);
        end
        chOn = false(1,nC); chOn(1) = true;
        % align hi-res pixel EDGES with base pixel EDGES -> exact registration
        set(hOv,'XData',[0.5+0.5*sy/W, sy+0.5-0.5*sy/W], ...
                'YData',[0.5+0.5*sx/H, sx+0.5-0.5*sx/H],'Visible','on');
        updateOverlay();
    end

    function updateOverlay()
        if isempty(ovNorm), return; end
        if ~ovVisible, set(hOv,'Visible','off'); setStatus('Overlay hidden'); return; end
        set(hOv,'Visible','on');
        lo = oCenter - oWidth/2; hi = oCenter + oWidth/2;
        if hi <= lo, hi = lo + 1e-3; end
        [H,W] = size(ovNorm{1});
        rgb = zeros(H,W,3,'single');
        for c = 1:numel(ovNorm)
            if chOn(c)
                ci = mod(c-1,size(chColors,1))+1;
                chn = min(max((ovNorm{c}-lo)./(hi-lo),0),1);
                for q = 1:3
                    if chColors(ci,q) > 0, rgb(:,:,q) = rgb(:,:,q) + chn*chColors(ci,q); end
                end
            end
        end
        rgb = min(rgb,1);
        set(hOv,'CData',rgb,'AlphaData',oAlpha*max(rgb,[],3));
        on = find(chOn);
        setStatus(sprintf('Overlay channels on: %s | opacity %.1f',mat2str(on),oAlpha));
    end

    function zoomView(f)
        xl = xlim(ax); yl = ylim(ax);
        cx = mean(xl); cy = mean(yl);
        nx = cx + (xl-cx)*f; ny = cy + (yl-cy)*f;
        nx = clampLim(nx,[0.5 sy+0.5]); ny = clampLim(ny,[0.5 sx+0.5]);
        if diff(nx) < 5 || diff(ny) < 5, return; end
        set(ax,'XLim',nx,'YLim',ny);
    end

    function panView(fx,fy)
        xl = xlim(ax); yl = ylim(ax);
        nx = clampLim(xl + fx*diff(xl),[0.5 sy+0.5]);
        ny = clampLim(yl + fy*diff(yl),[0.5 sx+0.5]);
        set(ax,'XLim',nx,'YLim',ny);
    end

    function L = clampLim(L,B)
        w = min(diff(L),diff(B));
        if L(1) < B(1), L = [B(1) B(1)+w]; end
        if L(2) > B(2), L = [B(2)-w B(2)]; end
    end

    function addROI()
        n = numel(masks)+1;
        setStatus(sprintf('Draw %s %d: click-drag around it, release to finish.',opt.Name,n));
        h = drawfreehand(ax,'Color','r','LineWidth',1.5);
        if isempty(h) || ~isvalid(h), return; end
        pos = h.Position; delete(h);
        if size(pos,1) < 3, setStatus('ROI too small - ignored.'); return; end
        m = poly2mask(pos(:,1),pos(:,2),sx,sy);
        if ~any(m(:)), setStatus('Empty ROI - ignored.'); return; end
        pi_ = 0; pf = 0;
        if hasParent, [pi_,pf] = roiParentAssignment(double(m),opt.ParentMask); end
        masks{n} = m; parentIdx(n,1) = pi_; parentFrac(n,1) = pf;
        if pi_ > 0, col = pcols(pi_,:); else, col = [1 1 1]; end
        B = bwboundaries(m);
        hLines{n} = plot(ax,B{1}(:,2),B{1}(:,1),'-','Color',col,'LineWidth',1.5);
        [r,c] = find(m);
        hTexts{n} = text(ax,mean(c),mean(r),num2str(n),'Color',col,'FontWeight','bold', ...
            'HorizontalAlignment','center');
        if hasParent && pi_ > 0
            setStatus(sprintf('%s %d is in %s %d (%.0f%% inside). Press d for next, Enter to finish.', ...
                opt.Name,n,opt.ParentName,pi_,100*pf));
        elseif hasParent
            setStatus(sprintf('%s %d is NOT inside any %s.',opt.Name,n,opt.ParentName));
        else
            setStatus(sprintf('%s %d added. Press d for next, Enter to finish.',opt.Name,n));
        end
        figure(fig);
    end

    function undoROI()
        n = numel(masks);
        if n == 0, return; end
        delete(hLines{n}); delete(hTexts{n});
        masks(n) = []; hLines(n) = []; hTexts(n) = [];
        parentIdx(n) = []; parentFrac(n) = [];
        setStatus(sprintf('Removed %s %d.',opt.Name,n));
    end

    function finish()
        setappdata(fig,'done',true);
    end

    function onClose(~,~)
        setappdata(fig,'done',true);
    end

    function setStatus(s)
        status = s;
        if isempty(opt.Title), t = sprintf('Draw %s ROIs',opt.Name); else, t = opt.Title; end
        if isvalid(ax), title(ax,{t,status},'Interpreter','none'); end
        drawnow limitrate;
    end

    function showHelp()
        fprintf(['\n=== %s ROI GUI keys ===\n' ...
            ' =/-   zoom in/out        arrows  pan       r  reset view\n' ...
            ' b/B   base bright+/-     n/N     base contrast+/-\n' ...
            ' o/O   overlay bright+/-  p/P     overlay contrast+/-\n' ...
            ' t/T   overlay opacity+/- 1,2,3.. toggle overlay channel   v  overlay on/off\n' ...
            ' d     draw next ROI      u  undo   Enter/f  FINISH\n'],opt.Name);
    end
end
