function [limits, handles] = interactivehistogram(x, y, v, m)
% [limits, handles] = interactivehistogram(x, y, v)
% [limits, handles] = interactivehistogram(x, y, v, m)
% [limits, handles] = interactivehistogram(x, y, v, xbins)
%
% input arguments v, m, xbins are the same as for hist(...), see help hist
%
% x, y, v are vectors of the same size
% display graph of x,y, and histogram of v. remove points from the x,y
% graph based on the thresholding of v
%
% limits = [min max] values of v as determined interactively with the cursors
% handles = struct with handles to the graphics

if ~exist('m','var') || isempty(m),
    m = 51;
end

% copy the inputs for later filtering
xin = x;
yin = y;

handles.fig = figure('position',[360 81 560 850]);
handles.axprof = subplot(2,1,1);
handles.profile = plot(x, y, '.-'); grid
handles.axhist = subplot(2,1,2); hist(v,m); grid

% create the cursor lines in the histogram plot, set their properties
xlimtmp = get(handles.axhist,'xlim');
ylimtmp = get(handles.axhist,'ylim');
handles.lineCursor(1) = line([0.7 0.3]*xlimtmp'*ones(2,1),ylimtmp);
handles.lineCursor(2) = line([0.3 0.7]*xlimtmp'*ones(2,1),ylimtmp);

set(handles.lineCursor(1),'linewidth',2,'color','b','ButtonDownFcn',@LineCursorCallback)
set(handles.lineCursor(2),'linewidth',2,'color','b','ButtonDownFcn',@LineCursorCallback)

% create variables to store return values
lev1 = [];
lev2 = [];

% create the callback functions:
    function LineCursorCallback(hh, evnt)
        set(hh,'Selected','on');
        set(handles.fig,...
            'WindowButtonMotionFcn', {@ButtonMotionCallback, hh}, ...
            'WindowButtonUpFcn', {@ButtonUpCallback, hh});

    end % LineCursorCallback

    function ButtonMotionCallback(hFig, evnt, hcursor)
        [x2,y2] = GetCurrentPt(handles.axhist);
        
        % limit how far user can slide cursors
        xlim = get(gca,'xlim');
        
        %         if x2 < min(v), x2 = min(v); end
        %         if x2 > max(v), x2 = max(v); end
        if x2 < min(xlim), x2 = min(xlim); end
        if x2 > max(xlim), x2 = max(xlim); end

        set(hcursor,'xdata',[x2; x2]);
        
        ApplyCursorThresholds;

    end % ButtonMotionCallback

    function ButtonUpCallback(hFig, evnt, hcursor)
        set(hFig, 'WindowButtonMotionFcn', '');
        set(hcursor, 'Selected','off');
        
        %PlotCursor(hcursor, hAxesVerCursor);

    end % ButtonUpCallback

    %%%%%%%%%%%%%%%%%%%%%%% what to do when the cursor moves %%%%%%%%%%%
    function ApplyCursorThresholds
        % get the cursor positions in the histogram plot, and use them
        % to threshold the profile data
        lev1 = get(handles.lineCursor(1),'xdata');
        lev2 = get(handles.lineCursor(2),'xdata');
        
        [x, y] = filterdata(...
            v > min([lev1(1) lev2(1)]) & ...
            v < max([lev1(1) lev2(1)]),...
            xin, yin);
        
        set(handles.profile,'xdata',x);
        set(handles.profile,'ydata',y);
        
    end % PlotCursor

uiwait(handles.fig);

% return the final limit values
limits = [lev1(1) lev2(1)];

end % main

%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CursorColor(hcursor, color)
set(hcursor,'color',color)
hPlot = get(hcursor,'userdata');
set(hPlot,'color',color);

end % CursorColor

function [x, y] = GetCurrentPt(hax)
    pt = get(hax, 'CurrentPoint');
    x = pt(1,1);
    y = pt(1,2);
end % GetCurrentPt
