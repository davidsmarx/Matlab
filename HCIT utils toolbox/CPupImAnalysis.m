 classdef CPupImAnalysis
    % S = CPupImAnalysis(pupilimage)
    %
    % grresponse = MaskEdges(S)
    % [psfedge, imedge, rline0, psffwhm] = EdgeDeconv(S, varargin)
    % Sm = CalcDimensions(S, varargin)
    
    
    properties
    
        U = CConstants;
        
        Im
        bMask
        MaskThresh
        camz
        pix  %= 6.5*1e-6; % um pixel size of pupil image (nbin = 1)
        x
        y
        X
        Y
        R
        T
        com
        pPad = 1.25; % factor
        sDims
        Sstruts

    end % properties
    
    
    methods
        function S = CPupImAnalysis(imfn, varargin)
    
            % options
            S.pix = CheckOption('pix', 6.5*S.U.UM, varargin{:});
            bDisplay = CheckOption('display', false, varargin{:});
            com = CheckOption('com', [], varargin{:}); % should be 1-offset for use as index
            CoordOffset = CheckOption('CoordOffset', [0 0], varargin{:}); % same units as pix
            Sstruts = CheckOption('Sstruts', [], varargin{:}); % if you want to analyze this image with struts already found with reference image
            
            
            % read image and header
            if isa(imfn, 'double')
                PupIm = imfn;
            else
                finfo = fitsinfo(PathTranslator(imfn));
                S.camz = FitsGetKeywordVal( finfo.PrimaryData.Keywords,'CAMZ');
                im = fitsread(PathTranslator(imfn),'image');
                PupIm = sqrt(im);
            end

            % center image and reduce padding
            [sResults, bMask] = AutoMetric(PupIm);
            S.MaskThresh = sResults.thresh;
            
            [nr, nc] = size(bMask);
            xx = 1:nc; yy = 1:nr; [XX, YY] = meshgrid(xx,yy);
            if isempty(com),
                com = [sum(XX(bMask)) sum(YY(bMask))]./sum(bMask(:));
            end
            S.com = round(com);
            
            %             figure_mxn(1,2)
            %             subplot(1,2,1), imageschcit(xx, yy, abs(PupIm))
            %             hold on
            %             plot(com(1), com(2), 'rx')
            %             hold off

            % limits for cropped image
            x1 = min(XX(bMask)); x2 = max(XX(bMask));
            y1 = min(YY(bMask)); y2 = max(YY(bMask));
            
            % width as even numbers
            wx = S.pPad*(x2-x1);
            wy = S.pPad*(y2-y1);
            wx = ceil(wx/2);
            wy = ceil(wy/2);

            % crop image and bMask
            iyuse = yy >= S.com(2) - wy & yy < S.com(2) + wy;
            ixuse = xx >= S.com(1) - wx & xx < S.com(1) + wx;
            S.Im = PupIm(iyuse, ixuse);
            S.bMask = bMask(iyuse, ixuse);

            %[S.x, S.y, S.X, S.Y, S.R, S.T] = CreateGrid(S.Im, S.pix);
            [x, y] = CreateGrid(S.Im, S.pix);
            S.x = x + CoordOffset(1);
            S.y = y + CoordOffset(2);
            [S.X, S. Y] = meshgrid(S.x, S.y);
            S.R = hypot(S.X, S.Y);
            S.T = atan2(S.Y, S.X);
            
            % dimensions
            S.sDims = S.CalcDimensions;
            
            if bDisplay,
                figure_mxn(1,2)
                subplot(1,2,1), imageschcit(S.x/S.U.UM,S.y/S.U.UM,abs(S.Im))
                title('|Image|')
                xlabel('Camera X (\mum)'), ylabel('Camera Y (\mum)')
                subplot(1,2,2), imageschcit(S.x/S.U.UM,S.y/S.U.UM,S.bMask)
                title(['Generated Binary Mask'])
                xlabel('Camera X (\mum)'), ylabel('Camera Y (\mum)')
            end
            
            if isempty(Sstruts),
                S.Sstruts = S.FindAllStruts(varargin{:});
            else
                S.Sstruts = Sstruts;
            end

            S.Sstruts.Spsf = S.DeconvAllStruts(varargin{:});
            %Spsf(istrut) = S.StrutDeconv(Sstrut(istrut), varargin{:});

            
        end % CPupImAnalysis
        
        function grresponse = MaskEdges(S)
            % don't think this is used anymore?
            
            [Fx, Fy] = gradient(double(S.bMask));
            Iedges = Fx.^2 + Fy.^2; % = 0, 0.25, or 0.5
            
            figure, imageschcit(S.x/S.U.UM, S.y/S.U.UM, Iedges)
            
            %bEdges = Iedges >= 0.5;
            
            %thmin = pi-0.1;       thmax = pi+0.1;
            thc = (-180+6:12:(180-6))*S.U.P;
            thmin = thc - 6*S.U.P;
            thmax = thc + 6*S.U.P;
            Nth = length(thc);
            
            nlines = 10;  % # of lines per sector
            nrline = 128; % # of pixels to resample edge responses
            grresponse = zeros(nrline, Nth);
            for ith = 1:Nth,
                th1 = thmin(ith);
                th2 = thmax(ith);
                
                %iuse = bEdges & (S.T > th1 & S.T <= th2) & (S.R > rmin & S.R < rmax);
                %meanR = mean(S.R( iuse ));

                rOD = max(S.R(:).*S.bMask(:).*(S.T(:) > th1 & S.T(:) < th2));
                rmin = rOD - 0.5*(S.U.MM);
                rmax = rOD + 0.5*(S.U.MM);
                
                thlines = linspace(th1, th2, nlines);
                rlines  = linspace(rmin,  rmax,  nrline);
                rline0  = rlines - mean(rlines); % centered around 0
                
                grlines = zeros(nrline,nlines);
                for jth = 1:nlines,
                    th = thlines(jth);
                    imline = interp2(S.X, S.Y, S.Im, rlines.*cos(th), rlines.*sin(th),'cubic',0);
                    grline = gradient(imline);
                    [~, xm] = findpeakfftinterp(rlines, grline, 3);
                    [~, ixm] = min(abs(xm - rOD)); % gradient peak closest to the binary edge
                    xm = xm(ixm);
                    
                    grlines(:, jth) = interp1(rlines - xm, grline, rline0, 'spline');
                end
                
                grresponse(:,ith) = mean(abs(grlines),2);
                
            end % for each theta band
            
            figure, plot(rline0/S.U.UM, grresponse), grid
                
            % choose an area to analyze
            
        end % MaskEdges

        function Sedge = EdgeDeconv(S, varargin)
                        
            % options
            bDebug = CheckOption('debug', false, varargin{:});
            thc = CheckOption('thcenter', 0.6*pi, varargin{:});
            thw = CheckOption('thwidth', 12*S.U.P, varargin{:});
            IDorOD = CheckOption('idorod', 'od', varargin{:});
            
            %
            [Fx, Fy] = gradient(double(S.bMask));
            Iedges = Fx.^2 + Fy.^2; % = 0, 0.25, or 0.5
            
            if bDebug,
                hfig = figure; imageschcit(S.x/S.U.UM, S.y/S.U.UM, Iedges)
            end

            % azimuthal location to use
            thmin = thc - thw/2;
            thmax = thc + thw/2;
            Nth = length(thc); % # number of theta bands

            % for resampling along radial across edge
            nlines = 10;  % # of lines per sector
            nrline = 256; % # of pixels to resample edge responses

            switch IDorOD,
                case 'od',
                    K = triu(ones(nrline)); % kernel for deconvolution with edge
                    fFindEdge = @max;
                case 'id'
                    K = tril(ones(nrline)); % kernel for deconvolution with edge
                    fFindEdge = @min;
                otherwise,
                    error(['unknown edge type ' IDorOD ]);
            end
            
            [psfedge, imedge] = deal(zeros(nrline, Nth));
            psffwhm = zeros(Nth,1);
            
            for ith = 1:Nth,
                th1 = thmin(ith);
                th2 = thmax(ith);
                
                %iuse = bEdges & (S.T > th1 & S.T <= th2) & (S.R > rmin & S.R < rmax);
                %meanR = mean(S.R( iuse ));

                rOD = fFindEdge(S.R(S.bMask & (S.T > th1 & S.T < th2)));
                rmin = rOD - 0.5*(S.U.MM);
                rmax = rOD + 0.5*(S.U.MM);
                
                thlines = linspace(th1, th2, nlines)';
                rlines  = linspace(rmin,  rmax,  nrline)';
                rline0  = rlines - mean(rlines); % centered around 0
                
                imline = zeros(nrline,nlines);
                psfe = zeros(nrline,nlines);                
                for jth = 1:nlines,
                    th = thlines(jth);
                    
                    if bDebug,
                        figure(hfig); hold on, plot(rlines.*cos(th)/S.U.UM, rlines.*sin(th)/S.U.UM, '-r'), hold off
                    end
                    
                    imline(:, jth) = interp2(S.X, S.Y, S.Im, rlines.*cos(th), rlines.*sin(th),'cubic',0);

                    %                     grline = gradient(imline);
                    %                     [~, xm] = findpeakfftinterp(rlines, grline, 3);
                    %                     [~, ixm] = min(abs(xm - rOD)); % gradient peak closest to the binary edge
                    %                     xm = xm(ixm);
                    
                    psfe(:, jth) = K \ abs(imline(:, jth));
                end
                
                imedge(:,ith) = mean(abs(imline),2);
                psfedge(:,ith) = mean(abs(psfe),2);
                
                % fwhm, need to trim to eliminate deconv effects at
                % first/last
                [rltmp, psftmp] = filterdata(rline0 > 0.5*rline0(1) & rline0 < 0.5*rline0(end), ...
                    rline0, psfedge(:,ith));
                try,
                    psffwhm(ith) = fwhm(rltmp(:), psftmp(:));
                catch
                    psffwhm(ith) = nan;
                end
                   
                
            end % for each theta band
            
            if bDebug, figure, plot(rline0/S.U.UM, psfedge), grid,
                xlabel('Radius (\mum)'), end

            % return struct
            Sedge = struct(...
                'psfedge', psfedge ...
                ,'imedge', imedge ...
                ,'rline0', rline0 ...
                ,'psffwhm', psffwhm ...
                );
                        
        end % EdgeDeconv

        function Sm = CalcDimensions(S, varargin)
            %             Sm.meanODradius = meanODradius;
            %             Sm.meanIDradius = meanIDradius;
            %             Sm.theta_OD = thODuse;
            %             Sm.radius_OD = rODuse;
            %             Sm.theta_ID = thIDuse;
            %             Sm.radius_ID = rIDuse;

            % options
            bDebug = CheckOption('debug', false, varargin{:});
            
            %
            [Fx, Fy] = gradient(double(S.bMask));
            Iedges = Fx.^2 + Fy.^2; % = 0, 0.25, or 0.5
            
            if bDebug,
                hfig = figure_mxn(1,2);
                haxim = subplot(1,2,1);
                imageschcit(S.x/S.U.UM, S.y/S.U.UM, Iedges)
            end
                    
            % alternative method for choosing ID and OD edge points
            [redges, xedges, yedges] = filterdata(Iedges >= 0.5, S.R, S.X, S.Y);
            [resort, xesort, yesort] = sortdata({redges, xedges, yedges});
            [cnt, rebin] = hist(resort, length(resort)/10);
            [~, locs] = findpeaks(cnt, rebin);
            deltabin = diff(rebin([1 2]));
            % locs are bin values of radius where there are a large # of edges
            
            % ID is first bin with peak
            [rIDuse, xIDuse, yIDuse] = filterdata(resort >= locs(1) - 2*deltabin & resort <= locs(1) + 2*deltabin, resort, xesort, yesort);
            thIDuse = atan2(yIDuse, xIDuse);
            
            % OD is last bin with peak
            %[rODuse, xODuse, yODuse] = filterdata(resort >= locs(end) - 2*deltabin & resort <= locs(end) + 2*deltabin, resort, xesort, yesort);
            [rODuse, xODuse, yODuse] = filterdata(resort >= locs(end) - 2*deltabin, resort, xesort, yesort);
            thODuse = atan2(yODuse, xODuse);
            
            if bDebug,
                %figure(hfig)
                axes(haxim);
                hold on
                plot(xODuse/S.U.UM, yODuse/S.U.UM, '.r')
                plot(xIDuse/S.U.UM, yIDuse/S.U.UM, '.c')
                hold off
                subplot(1,2,2)
                [haxrad, hbb, hcc]  = plotyy(thIDuse/pi, rIDuse/S.U.UM, thODuse/pi, rODuse/S.U.UM);
                grid on
                set(hbb,'Marker','.','LineStyle','none')
                set(hcc,'Marker','.','LineStyle','none')                
            end
            
            meanODradius = mean(rODuse);
            meanIDradius = mean(rIDuse);

            % result struct
            Sm.meanODradius = meanODradius;
            Sm.meanIDradius = meanIDradius;
            Sm.theta_OD = thODuse;
            Sm.radius_OD = rODuse;
            Sm.theta_ID = thIDuse;
            Sm.radius_ID = rIDuse;
            
        end % CalcDimensions
        
        function StrutObjects = FindAllStruts(S, varargin)
            % use bwboundaries to locate the struts
            % return strutxy1 and strutxy2 for use with FindStrut()
            % call FindStrut() for each strut to find edges of each strut
            % call StrutDeconv for each strut, edge gradient analysis
            
            rpad = CheckOption('radiuspad', 10*S.pix, varargin{:}); % um
            bDebug = CheckOption('debug', false, varargin{:});
            useObjects = CheckOption('useobjects', [], varargin{:});
            strutlengthfactor = CheckOption('strutlengthfactor', 1, varargin{:});
            %S.sDims = S.CalcDimensions;
            
            % first manipulate mask so that struts are 1, everywhere else 0
            bwmask = ~S.bMask;
            bwmask(S.R > S.sDims.meanODradius - rpad) = 0;
            bwmask(S.R < S.sDims.meanIDradius + rpad) = 0;
            
            % find objects
            [B,L,N,A] = bwboundaries(bwmask, 'noholes');
            % eliminate isolated pixels, etc
            for il = 1:N,
                objarea(il) = bwarea(L==il);
                if objarea(il) < 100,
                    %S.bMask(L==il) = true;
                    bwmask(L==il) = 0;
                end
            end            
            [B,L,N,A] = bwboundaries(bwmask, 'noholes');
        
            % use only requested objects
            if isempty(useObjects),
                useObjects = 1:N;
            else
                N = length(useObjects);
            end
            
            %             if N ~= 6,
            %                 figure, imageschcit(L), colorbar
            %                 error('wrong number of struts');
            %             end
        
            [strutxy1, strutxy2] = deal(zeros(N,2));
            for iobj = 1:N,
                istrut = useObjects(iobj);
                xstr = S.X(L==istrut);
                ystr = S.Y(L==istrut);
                rstr = S.R(L==istrut);
                
                %r1 = S.sDims.meanODradius - rpad;% - 1*S.U.MM;
                %                 [rtmp, xtmp, ytmp] = filterdata(rstr > r1 - rpad & rstr < r1 + rpad,...
                %                     rstr, xstr, ystr);
                [rtmp, xtmp, ytmp] = filterdata(rstr > max(rstr) - 0.1*range(rstr), ...
                    rstr, xstr, ystr);                
                strutxy1(iobj,:) = [mean(xtmp) mean(ytmp)];

                %                 r2 = S.sDims.meanIDradius + 1*S.U.MM;
                %                 [rtmp, xtmp, ytmp] = filterdata(rstr > r2 - rpad & rstr < r2 + rpad,...
                %                     rstr, xstr, ystr);
                [rtmp, xtmp, ytmp] = filterdata(rstr < min(rstr) + 0.1*range(rstr), ...
                    rstr, xstr, ystr);
                strutxy2(iobj,:) = [mean(xtmp) mean(ytmp)];
                                                       
            end % for each strut
            
            if bDebug,
                figure, imageschcit(S.x/S.U.MM, S.y/S.U.MM, L), colorbar
                for ii = 1:N,
                   text(strutxy1(ii,1)/S.U.MM, strutxy1(ii,2)/S.U.MM, ['strut #' num2str(ii)], 'FontSize', 10, 'Color', 'w');
                end
                hold on
                plot(strutxy1(:,1)/S.U.MM, strutxy1(:,2)/S.U.MM, '+w', ...
                    strutxy2(:,1)/S.U.MM, strutxy2(:,2)/S.U.MM, 'ow')
                hold off
            end

            % call FindStrut and Deconv for each strut
            if bDebug, hfigDebug = figure; imageschcit(S.x/S.U.UM, S.y/S.U.UM, S.Im); else, hfigDebug = []; end
            for istrut = 1:N,
                % how much of the strut to use?
                xyc = mean([strutxy1(istrut,:); strutxy2(istrut,:)]);
                xy1tmp = xyc + strutlengthfactor*(strutxy1(istrut,:)-xyc);
                xy2tmp = xyc + strutlengthfactor*(strutxy2(istrut,:)-xyc);
            
                Sstrut(istrut) = S.FindStrut(...
                    xy1tmp, xy2tmp ... strutxy1(istrut,:), strutxy2(istrut,:) ...
                    , 'hfigDebug', hfigDebug, varargin{:} ...
                    );
                
                %Spsf(istrut) = S.StrutDeconv(Sstrut(istrut), varargin{:});
                
            end % for eah strut
            
            
            StrutObjects = struct(...
                'N', N ...
                ,'strutxy1', strutxy1 ...
                ,'strutxy2', strutxy2 ...
                ,'Sstrut', Sstrut ...
                ...,'Spsf', Spsf ...
                ,'radiuspad', rpad ...
                );
            
        end % FindAllStruts
        
        function Spsf = DeconvAllStruts(S, varargin)

            for istrut = 1:S.Sstruts.N,
                Spsf(istrut) = S.StrutDeconv(S.Sstruts.Sstrut(istrut), varargin{:});
            end

        end % DeconvAllStruts
        
        function Sout = FindStrut(S, xy0, xy1, varargin)
            % find left and right edges of strut, given two points in the
            % strut gap
            %
            % xy0, xy1 are in (m) = (S.pix*pixels)
            %
            % options:
            %   'widthd' = width (m) of line across strut to find edges,
            %              default = 750um
            %   'nlines' = # of lines across strut to use to find edges,
            %              spaced between xy0 and xy1, default = 10
            %   'debug'
            
            widthd = CheckOption('widthd', 750*S.U.UM, varargin{:});
            nlines = CheckOption('nlines', 10, varargin{:});
            bDebug = CheckOption('debug', false, varargin{:});
            hfigDebug = CheckOption('hfigDebug', [], varargin{:});
            
            % lines perpendicular to initial approximate
            dxy  = diff([xy0; xy1]); % slope m  = dxy(2)/dxy(1)
            dxyp = [1 -1].*fliplr(dxy);     % perpendicular slope mp = -1/m
            % for some distance d from xy0 along perpendicular line with
            % slope dxyp, see notebook p.84
            d = linspace(-0.5*widthd,0.5*widthd)'; % (m)
            ddx = dxyp(1)./sqrt(dxyp*dxyp'); % dx = ddx * d
            ddy = dxyp(2)./sqrt(dxyp*dxyp');    % dy = ddy * d
            
            % now create perpendicular lines going down the strut from xy0
            % to xy1
            fStrutLine = @(x, m) m*x + xy0(2) - m*xy0(1);
            
            if bDebug,
                hfig = figure(101);
                figure_mxn(hfig,2,1);
                hax(1) = subplot(2,1,1);
                imageschcit(S.x/S.U.UM, S.y/S.U.UM, S.Im)
                title('FindStrut')
                hold on
                plot([xy0(1) xy1(1)]/S.U.UM, [xy0(2) xy1(2)]/S.U.UM, '-r')
                
                hax(2) = subplot(2,1,2);
            end
            
            % collect edge points
            [xlt, ylt, xrt, yrt] = deal(zeros(nlines,1));
            for ii = 1:nlines,
                xs = xy0(1) + ( (ii-1)/(nlines-1) )*dxy(1);
                ys = fStrutLine(xs, dxy(2)./dxy(1));
                
                xp = xs + ddx.*d; % xy0(1) + dx;
                yp = ys + ddy.*d; % xy0(2) + dy;
               
                % interpolate image along line, gradient, find edges
                imline = interp2(S.X, S.Y, S.Im, xp, yp);
                grimline = gradient(imline);
                % left edge is max negative gradient
                [glt, dlt] = findpeakn(d(d<0),-grimline(d<0));
                [grt, drt] = findpeakn(d(d>0), grimline(d>0));

                xlt(ii) = xs + ddx.*dlt;
                ylt(ii) = ys + ddy.*dlt;
                xrt(ii) = xs + ddx.*drt;
                yrt(ii) = ys + ddy.*drt;

                if bDebug,
                    axes(hax(1));
                    plot(xp/S.U.UM, yp/S.U.UM, '-b', ...
                        xlt(ii)/S.U.UM, ylt(ii)/S.U.UM, '.r', ...
                        xrt(ii)/S.U.UM, yrt(ii)/S.U.UM, '.g')
                    axes(hax(2));
                    [haxyy, h1 h2] = plotyy(d/S.U.UM, imline, d/S.U.UM, grimline);
                    xlabel('d (\mum)')
                    grid on
                    hold(haxyy(2),'on')
                    plot(haxyy(2), dlt*[1 1]/S.U.UM, get(haxyy(2),'ylim'), '--r', ...
                        drt*[1 1]/S.U.UM, get(haxyy(2),'ylim'), '--r')
                    grid on
                end
                
            end % 
            
            % linear fit to detected edges
            [p, s, mu] = polyfit(xlt, ylt, 1);
            yltfit = polyval(p, xlt, [], mu);
            [p, s, mu] = polyfit(xrt, yrt, 1);
            yrtfit = polyval(p, xrt, [], mu);
            
            if bDebug,
                if isempty(hfigDebug),
                    hfigDebug = figure;
                    imageschcit(S.x/S.U.UM, S.y/S.U.UM, S.Im)
                else
                    figure(hfigDebug),                    
                end
                hold on
                plot(xlt/S.U.UM, ylt/S.U.UM, '.r', xrt/S.U.UM, yrt/S.U.UM, '.r', ...
                    xlt/S.U.UM, yltfit/S.U.UM, '-r', xrt/S.U.UM, yrtfit/S.U.UM, '-r')
            end
            
            % angle between left and right edges
            thlt = atan2(diff(yltfit([1 end])), diff(xlt([1 end])));
            thrt = atan2(diff(yrtfit([1 end])), diff(xrt([1 end])));
            
            % mean distance between left and right edges
            StrutWidth = sqrt( (xlt-xrt).^2 + (ylt-yrt).^2 );
            meanStrutWidth = mean(StrutWidth);
            
            % center of strut
            xStrutCenter = mean([xlt xrt],2);
            yStrutCenter = mean([ylt yrt],2);
            
            Sout = struct(...
                'xlt',xlt ...
                ,'ylt',ylt ...
                ,'xrt',xrt ...
                ,'yrt',yrt ...
                ,'yltfit', yltfit ...
                ,'yrtfit', yrtfit ...
                ,'thlt', thlt ...
                ,'thrt', thrt ...
                ,'meanStrutWidth', meanStrutWidth ...
                ,'xStrutCenter', xStrutCenter ...
                ,'yStrutCenter', yStrutCenter ...
                );
            
        end % FindStrut
    
        function Spsf = StrutDeconv(S, Sstrut, varargin)
            % use Kernel with two edges to deconv strut image to get 1-d
            % psf
            % Sstrut is struct returned by method FindStruct
        
            widthd = CheckOption('widthd', 750*S.U.UM, varargin{:}); % length of cross-section line; not used except to create Ns
            bOuterEdges = CheckOption('outeredges', false, varargin{:}); % look for outside pair of edges; for edge mask analysis
            
            
            ds = S.pix/8; % sample spacing for resampling and deconv kernel
            %Ls = 4*Sstrut.meanStrutWidth; % not used except to create Ns

            
            bDebug = CheckOption('debug', false, varargin{:});
            
            Ns = ceil(widthd/ds); % # of samples across strut
            d  = CreateGrid(Ns, ds);
            
            Nlines = length(Sstrut.xStrutCenter);
            
            % lines perpendicular to line along center of strut, the lines
            % calculated in FindStrut were an estimate, these should be
            % more accurate      
            % end points:
            xy0 = [Sstrut.xStrutCenter(1)   Sstrut.yStrutCenter(1)];
            xy1 = [Sstrut.xStrutCenter(end) Sstrut.yStrutCenter(end)];
            
            dxy  = diff([xy0; xy1]); % slope m  = dxy(2)/dxy(1)
            dxyp = [1 -1].*fliplr(dxy);     % perpendicular slope mp = -1/m
            % for some distance d from xy0 along perpendicular line with
            % slope dxyp, see notebook p.84
            ddx = dxyp(1)./sqrt(dxyp*dxyp'); % dx = ddx * d
            ddy = dxyp(2)./sqrt(dxyp*dxyp'); % dy = ddy * d

            % create Kernel of two edges
            % first estimate of strut width
            nWidth = round(Sstrut.meanStrutWidth/ds);
            KK = toeplitz([ones(1,Ns-nWidth) zeros(1,nWidth) ones(1,Ns-nWidth)]);
            K = KK(1:Ns, Ns-nWidth+1:end); % Ns x Ns
            
            if bDebug,
                hfig = figure(102);
                figure_mxn(hfig,2,1)
                hax = [subplot(2,1,1) subplot(2,1,2)];
                axes(hax(1));
                imageschcit(S.x/S.U.UM, S.y/S.U.UM, S.Im)
                title('StrutDeconv')

            end
            
            Kc = [];
            rhs = [];
            [grAs, As] = deal(zeros(Ns, Nlines));
            for ii = 1:Nlines,
                xs = Sstrut.xStrutCenter(ii) + ddx*d;
                ys = Sstrut.yStrutCenter(ii) + ddy*d;
                
                % image amplitude along cut
                As(:,ii) = interp2(S.X, S.Y, sqrt(S.Im), xs, ys);

                % deconv (decided not to use this for anything)
                psfii = K \ As(:,ii);
                
                % psf from gradient
                grAs(:,ii) = gradient(As(:,ii));
                % location and value of gradient peaks
                [grMaxLt(ii), xgrMaxLt(ii), ygrMaxLt(ii)] = grmin(grAs(:,ii),... % bright-dark edge is negative gradient
                    d<0 & d > -1.5*Sstrut.meanStrutWidth, xs, ys); % used to be & d > -1.5*Sstrut.meanStrutWidth
                [grMaxRt(ii), xgrMaxRt(ii), ygrMaxRt(ii)] = grmax(grAs(:,ii),... % dark-bright edge is positive
                    d>0 & d <  1.5*Sstrut.meanStrutWidth, xs, ys);
                
                if bOuterEdges,
                    
                    [grMaxLt_out(ii), xgrMaxLt_out(ii), ygrMaxLt_out(ii)] = grmax( grAs(:,ii),...
                        d<0, xs, ys);
                    [grMaxRt_out(ii), xgrMaxRt_out(ii), ygrMaxRt_out(ii)] = grmin( grAs(:,ii),...
                        d>0, xs, ys);
                else
                    
            
                end
                
                if bDebug,
                    axes(hax(1));
                    hold on
                    plot(xs/S.U.UM, ys/S.U.UM, '-r')
                    hold off
                    
                    axes(hax(2));
                    plot(d/S.U.UM, As(:,ii) ...
                        ,d/S.U.UM, grAs(:,ii) ...
                        );    %                         d/S.U.UM, K(floor((Ns-nWidth)/2),:), ...
                    %                         d/S.U.UM, psfii ...
                        
                    grid on
                end
            
                Kc = [Kc; K];
                rhs = [rhs; As(:,ii)];
            end % for each cross-strut line
            
            % deconvolve, doesn't work too well
            %psfpsf = Kc \ rhs;
            
            % mean edge gradient
            meangrAs = mean(grAs,2);
            %             % peak gradient values for left and right edges
            %             [peakgrLt, dpeakgrLt] = grmax(-meangrAs, d<0 & d>(-1.5*Sstrut.meanStrutWidth), d);
            %             [peakgrRt, dpeakgrRt] = grmax( meangrAs, d>0 & d<( 1.5*Sstrut.meanStrutWidth), d);
            
            [fwhmlt, ~, xelt, dhalflt] = fwhm(d(d < 0), -meangrAs(d < 0));
            [fwhmrt, ~, xert, dhalfrt] = fwhm(d(d > 0),  meangrAs(d > 0));
            
            % if there are outside edges to the left and right bright
            % areas, for example the pupil image quality edge test mask has
            % dark-bright edges outside the strut edges
            if bOuterEdges,
                [fwhmlt_out, ~, xelt_out, dhalflt_out] = fwhm(d(d < 0),  meangrAs(d < 0));
                [fwhmrt_out, ~, xert_out, dhalfrt_out] = fwhm(d(d > 0), -meangrAs(d > 0));
            else
                [fwhmlt_out, xelt_out, dhalflt_out, fwhmrt_out, xert_out, dhalfrt_out] = deal([]);
            
            end
            
            strutwidth = abs(xert - xelt);

            Spsf = struct(...
                'Sstrut', Sstrut ...
                ,'d', d ...
                ,'ddx', ddx ...
                ,'ddy', ddy ...
                ,'As', As ...
                ,'grAs', grAs ...
                ,'meangrAs', meangrAs ...
                ...%                 ,'peakgrLt', peakgrLt ...
                ...%                 ,'peakgrRt', peakgrRt ...
                ,'grMaxLt', grMaxLt ...
                ,'xgrMaxLt', xgrMaxLt ...
                ,'ygrMaxLt', ygrMaxLt ...
                ,'grMaxRt', grMaxRt ...
                ,'xgrMaxRt', xgrMaxRt ...
                ,'ygrMaxRt', ygrMaxRt ...
                ,'fwhmlt', fwhmlt ...
                ,'fwhmrt', fwhmrt ...
                ,'dhalflt', dhalflt ...
                ,'dhalfrt', dhalfrt ...
                ,'strutwidth', strutwidth ...
                ,'grMaxRt_out', grMaxLt_out...
                ,'grMaxLt_out', grMaxRt_out...
                ,'fwhmlt_out', fwhmlt_out ...
                ,'fwhmrt_out', fwhmrt_out ...
                ,'dhalflt_out', dhalflt_out ...
                ,'dhalfrt_out', dhalfrt_out ...
                );
            
        end % StrutDeconv
        
        function [hfig, hax] = DisplayImage(S, varargin)
            
            hfig = figure;
            imageschcit(S.x/S.pix, S.y/S.pix, S.Im)
            xlabel('Camera X (pixel)')
            ylabel('Camera Y (pixel)')
            
        end % DisplayImage
        
        function [hfig, hax] = DisplayMask(S, varargin)
        end % DisplayMask
        
        function [hfig, hax] = PlotEdgeResponse(S, Spsf, varargin)
            % [hfig, hax] = PlotEdgeResponse(S, Spsf, varargin)
                        
            hax = CheckOption('hax', [], varargin{:});
            xlim = CheckOption('xlim', [], varargin{:});
            stitle = CheckOption('title', '', varargin{:});
            
            % [fwhmlt, ~, xelt, dhalflt] = fwhm(Spsf.d(Spsf.d < 0), -Spsf.meangrAs(Spsf.d < 0));
            % [fwhmrt, ~, xert, dhalfrt] = fwhm(Spsf.d(Spsf.d > 0),  Spsf.meangrAs(Spsf.d > 0));
            % strutwidth = abs(xert - xelt);
            
            if isempty(hax),
                hfig = figure;
                hax = gca;
            else
                axes(hax);
            end
            
            plot(Spsf.d/S.U.UM, [Spsf.grAs]), grid
            ylim = get(gca,'ylim');
            hold on
            hpl = plot(Spsf.d/S.U.UM, Spsf.meangrAs, '.r');
            hpl.LineWidth = 1;
            plot(Spsf.dhalflt(1)*[1 1]/S.U.UM, ylim, '--r', Spsf.dhalflt(2)*[1 1]/S.U.UM, ylim, '--r')
            plot(Spsf.dhalfrt(1)*[1 1]/S.U.UM, ylim, '--r', Spsf.dhalfrt(2)*[1 1]/S.U.UM, ylim, '--r')
            text(Spsf.dhalflt(2)/S.U.UM, 0.75*ylim(1), [' ' num2str(Spsf.fwhmlt/S.U.UM,'%.1f') '\mum'],'fontsize',12,'color','r')
            text(Spsf.dhalfrt(2)/S.U.UM, 0.75*ylim(2), [' ' num2str(Spsf.fwhmrt/S.U.UM,'%.1f') '\mum'],'fontsize',12,'color','r')

            if ~any([isempty(Spsf.dhalflt_out) isempty(Spsf.dhalfrt_out)])
                plot(Spsf.dhalflt_out(1)*[1 1]/S.U.UM, ylim, '--r', Spsf.dhalflt_out(2)*[1 1]/S.U.UM, ylim, '--r')
                plot(Spsf.dhalfrt_out(1)*[1 1]/S.U.UM, ylim, '--r', Spsf.dhalfrt_out(2)*[1 1]/S.U.UM, ylim, '--r')
                text(Spsf.dhalflt_out(2)/S.U.UM, 0.75*ylim(2), [' ' num2str(Spsf.fwhmlt_out/S.U.UM,'%.1f') '\mum'],'fontsize',12,'color','r')
                text(Spsf.dhalfrt_out(2)/S.U.UM, 0.75*ylim(1), [' ' num2str(Spsf.fwhmrt_out/S.U.UM,'%.1f') '\mum'],'fontsize',12,'color','r')                
            end
            
            xlabel('Cross Section (\mum)'), ylabel('Gradient Intensity')
            if ~isempty(xlim), set(gca,'xlim',xlim/S.U.UM), end
            %set(gca,'ylim',ylim)
            
            title([stitle ', ' num2str(Spsf.strutwidth/S.U.UM,'%.1f') '\mum'])
            
        end % PlotEdgeResponse
        
        function [hfig, hax, han] = DisplayEdgeResponseAllStruts(S, varargin)
            % [hfig, hax, han] = DisplayEdgeResponseAllStruts(S, varargin)                      
            
            titlestr = CheckOption('title', [], varargin{:});
            xlim = CheckOption('xlim', [], varargin{:});
            
            hfig = figure_mxn(ceil(S.Sstruts.N/2), 2);
            for istrut = 1:S.Sstruts.N,
                hax = subplot(ceil(S.Sstruts.N/2), 2, istrut);
                S.PlotEdgeResponse(S.Sstruts.Spsf(istrut),'xlim', xlim,...
                    'hax', hax, 'title', ['strut #' num2str(istrut)])
        
            end % for istrut
    
            % annotation title at top
            han = annotation('textbox', [0.5 0.8 0.2 0.2],'String', titlestr ,...
                'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 16, 'Color', 'r', ...
                'FontWeight', 'bold');
            set(han,'HorizontalAlignment','center')
            ppp = get(han,'Position');
            set(han,'Position',[0.5 - 0.5*ppp(3) ppp(2:end)])
            
        end % DisplayEdgeResponseAllStruts
        
    end % methods
    
    
end % classdef

function [peakgr, varargout] = grmax(gr, iuse, varargin)

    cArgout = filterdata(iuse, gr, varargin{:});
    [peakgr, imax] = max(cArgout{1});
    %dpeak = dtmp(imax);
    for ii = 1:length(varargin)
        varargout{ii} = varargin{ii}(imax);
    end
    
end % grmax

function [peakgr, varargout] = grmin(gr, iuse, varargin)

    cArgout = filterdata(iuse, gr, varargin{:});
    [peakgr, imax] = min(cArgout{1});
    %dpeak = dtmp(imax);
    for ii = 1:length(varargin)
        varargout{ii} = varargin{ii}(imax);
    end
    
end % grmin