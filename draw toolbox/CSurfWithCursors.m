function S = CSurfWithCursors(x, y, z)
% S = CSurfWithCursors(x, y, z)

% % regrid with griddata
% xr = linspace(min(x(:)),max(x(:)),51);
% yr = linspace(min(y(:)),max(y(:)),51)';
% [XI, YI, ZI] = griddata(x, y, z, xr, yr,'linear');
[XI, YI, ZI] = deal(x, y, z);
xr = XI(1,:);
yr = YI(:,1);

% set up figure and axes size and position
screensize = get(0,'ScreenSize');
screenwidth = screensize(3);
screenheight = screensize(4);
figheight = screenheight/2;
figwidth = 3*figheight;
hFig = figure('position',[0 screenheight-figheight-75 figwidth figheight]);
SetFigureMenuOptions(hFig);

% set up three axes for the image surf and two cursor graphs
axeshscale = figwidth/figheight;
hAxesImage = axes('position',[0.1 0.15 0.8/axeshscale 0.8]);
xtmp = 0.1 + 0.8/axeshscale + 0.1;
hAxesHorCursor = axes('position',[xtmp 0.6 0.9-xtmp 0.3]);
hAxesVerCursor = axes('position',[xtmp 0.15 0.9-xtmp 0.3]);

% display the image in its axes
hSurf = surf(hAxesImage,xr,yr,ZI,...
    'linestyle','none',...
    'FaceColor','interp');
%hSurf = mesh(hAxesImage,xr,yr,ZI);

% set(hSurf,'edgecolor','none');
% set(hSurf,'FaceColor','interp');
% axis(hAxesImage,'image');
% set(hSurf,'alphaData',alphamap);
grid(hAxesImage,'on');
hColorbar = colorbar('peer',hAxesImage);
colorbartitle(hColorbar,'Z (\mum)')

% the cursors are planes, add patch objects to current axes
axes(hAxesImage);
xlimtmp = get(hAxesImage,'xlim');
ylimtmp = get(hAxesImage,'ylim');
zlimtmp = get(hAxesImage,'zlim');
pprops = {'ButtonDownFcn',@LineCursorCallback,'FaceAlpha',0.1};
hPatchVerCursor(1) = patch([0.7 0.3]*xlimtmp'*ones(4,1), ylimtmp([1 2 2 1])', zlimtmp([1 1 2 2])','b',pprops{:});
hPatchVerCursor(2) = patch([0.3 0.7]*xlimtmp'*ones(4,1), ylimtmp([1 2 2 1])', zlimtmp([1 1 2 2])','b',pprops{:});
hPatchHorCursor(1) = patch(xlimtmp([1 2 2 1])', [0.7 0.3]*ylimtmp'*ones(4,1), zlimtmp([1 1 2 2])','r',pprops{:});
hPatchHorCursor(2) = patch(xlimtmp([1 2 2 1])', [0.3 0.7]*ylimtmp'*ones(4,1), zlimtmp([1 1 2 2])','r',pprops{:});

% initialize, plot each cursor
% create a lineseries object in the appropriate axes for each cursor
% link the cursor plot to the cursor through its userdata
set(hPatchVerCursor(1),'userdata',InitialPlotCursor(hPatchVerCursor(1), hAxesVerCursor));
set(hPatchVerCursor(2),'userdata',InitialPlotCursor(hPatchVerCursor(2), hAxesVerCursor));
set(hPatchHorCursor(1),'userdata',InitialPlotCursor(hPatchHorCursor(1), hAxesHorCursor));
set(hPatchHorCursor(2),'userdata',InitialPlotCursor(hPatchHorCursor(2), hAxesHorCursor));

set(hAxesImage,'xlimmode','manual','ylimmode','manual','zlimmode','manual');

    %%%%%%%%%%%%%%%%%%%%%%% cursor plotting functions %%%%%%%%%%%
    function hPlot = InitialPlotCursor(hcursor, hax)
        % initial creation of the lineseries object for hcursor in hax
        
        axes(hax)
        hold on;
        
        xx = get(hcursor,'xdata');
        yy = get(hcursor,'ydata');

        % is this a horizontal or vertical cursor?
        if any(hcursor == hPatchHorCursor),
            % horizontal cursor, find the index of the y value
            [dtmp, iy] = min(abs(yr - yy(1)));
            zz = ZI(iy, :);
            % plot to the cursor graph
            hPlot = plot(hax,xr,zz);
            set(hPlot,'color',get(hcursor,'FaceColor'));
        else
            % vertical cursor
            [dtmp, ix] = min(abs(xr - xx(1)));
            zz = ZI(:, ix);
            hPlot = plot(hax,yr,zz);
            set(hPlot,'color',get(hcursor,'FaceColor'));
        end

        hold off;
        grid(hax,'on');
        
    end % InitialPlotCursor

    %%%%%%%%%%%%%%%%%%%%%%% cursor plotting functions %%%%%%%%%%%
    function PlotCursor(hcursor)
        % get the data for hcursor
        % modify the plot data for the plot associated with this cursor
        
        %axes(hax)
        hPlot = get(hcursor,'userdata');
        
        xx = get(hcursor,'xdata');
        yy = get(hcursor,'ydata');

        % is this a horizontal or vertical cursor?
        if any(hcursor == hPatchVerCursor),
            % vertical cursor
            [dtmp, ix] = min(abs(xr - xx(1)));
            zz = ZI(:, ix);
            xx = yr;
            yy = zz;
        else
            % horizontal cursor, find the index of the y value
            [dtmp, iy] = min(abs(yr - yy(1)));
            zz = ZI(iy, :);
            xx = xr;
            yy = zz;
        end
        
        set(hPlot,'xdata',xx);
        set(hPlot,'ydata',yy);
        
    end % PlotCursor

%%%%%%%%%%%%%%%%%%%%%%%% CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%%%%
InitialCursorPosition = [];
ButtonDownMousePosition = [];
    function LineCursorCallback(hh, evnt)
        set(hh,'Selected','on');
        
        % record initial cursor location and mouse location
        InitialCursorPosition = [get(hh,'xdata') get(hh,'ydata')];
        ButtonDownMousePosition = GetCurrentPt(hAxesImage);
        
        set(hFig,...
            'WindowButtonMotionFcn', {@ButtonMotionCallback, hh}, ...
            'WindowButtonUpFcn', {@ButtonUpCallback, hh});

    end % LineCursorCallback

    function ButtonMotionCallback(hFig, evnt, hcursor)
        xypt = GetCurrentPt(hAxesImage);
                
        % check that current point is within the axes limits
        xlim = get(hAxesImage,'xlim');
        ylim = get(hAxesImage,'ylim');
        if xypt(1) >= xlim(1) & xypt(1) <= xlim(2) & xypt(2) >= ylim(1) & xypt(2) <= ylim(2),
            % is this a vertical or horizontal cursor
            if any(hcursor == hPatchVerCursor)
                xnew = InitialCursorPosition(1,1) + (xypt(1) - ButtonDownMousePosition(1));
                set(hcursor,'xdata',xnew*ones(4,1));
            else
                ynew = InitialCursorPosition(1,2) + (xypt(2) - ButtonDownMousePosition(2));
                set(hcursor,'ydata',ynew*ones(4,1))
            end

            PlotCursor(hcursor);

        end % if current point inside axes
        
    end % ButtonMotionCallback

    function ButtonUpCallback(hFig, evnt, hcursor)
        set(hFig, 'WindowButtonMotionFcn', '');
        set(hcursor, 'Selected','off');
        
        %PlotCursor(hcursor, hAxesVerCursor);

    end % ButtonUpCallback

S.hFig = hFig;
S.hAxesImage = hAxesImage;
S.hSurf = hSurf;
S.hColorbar = hColorbar;
S.hAxesHorCursor = hAxesHorCursor;
S.hAxesVerCursor = hAxesVerCursor;
S.hPatchHorCursor = hPatchHorCursor;
S.hPatchVerCursor = hPatchVerCursor;
S.CursorColor = @CursorColor;

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CursorColor(hcursor, color)
set(hcursor,'FaceColor',color)
hPlot = get(hcursor,'userdata');
set(hPlot,'color',color);

end % CursorColor


function xy = GetCurrentPt(hax)
    pt = get(hax, 'CurrentPoint');
    x = pt(1,1);
    y = pt(1,2);
    xy = [x y];
end % GetCurrentPt