classdef CDataSet4D < handle
    % S = CDataSet4D
    % S = CDataSet4D([nrows ncols]) creates and array of empty CDataSet4D's
    % S = CDataSet4D(x, y, z, zname, w, wname)
    % S = CDataSet4D(x, y, z, zname)
    % methods:
    %   copy
    %   plotz
    %   removetilt
    %   is4D
    %   SetHolevalNaN
    %   ApplyUnits
    %   FilterData
    %   cat
    %   OutlierFilter
    %   Resample
    %   filtfilt

    properties
        x
        y
        z
        zname = '';
        w
        wname = '';
        tag = '';
        hParent;
        % holeval = NaN;
        
    end % properties
    
    events
        %eventCheckTiltCorrection
        
    end % events
    
    methods
        function S = CDataSet4D(varargin)
            switch nargin,
                case 0,
                    % S = CDataSet
                    
                case 1,
                    nm = varargin{1};
                    % doesn't work because repmat copies the handle:
                    %S = repmat(CDataSet4D, nm(1), nm(2));
                    for ii = 1:nm(1)*nm(2),
                        S(ii) = CDataSet4D;
                    end
                    S = reshape(S, nm(1), nm(2));

                case 4,
                    % S = CDataSet(x, y, z, zname)
                    [S.x, S.y, S.z, S.zname] = deal(varargin{:});

                case 6,
                    % S = CDataSet(x, y, z, zname, w, wname)
                    [S.x, S.y, S.z, S.zname, S.w, S.wname] = deal(varargin{:});

                otherwise,
                    error('usage');
            end % switch nargin
               
        end % CDataSet
        
        function Scopy = copy(S)
            Scopy = CDataSet4D(S.x, S.y, S.z, S.zname, S.w, S.wname);
            Scopy.tag = S.tag;
%             Scopy.tiltcorrection = S.tiltcorrection;
            
        end % copy
               
        function h = plotz(S, Units)
            % Units for graph

            if exist('Units','var')
                switch class(Units)
                    case 'struct'
                        [xunits, yunits, zunits] = deal(Units.xunits, Units.yunits, Units.zunits);
                    case 'double'
                        [xunits, yunits, zunits] = deal(Units(1),Units(2),Units(3));
                    case 'cell'
                        [xunits, yunits, zunits] = deal(Units{:});
                    otherwise,
                        error('units');
                end % switch
            else
                [xunits, yunits, zunits] = deal(1);
            end
            
            if length(S) > 1,
                stmp = S.cat;
            else
                stmp = S;
            end
            
            if is4D(stmp),
                h = plot3c(stmp.x/xunits, stmp.y/yunits, stmp.z/zunits, stmp.w, '.');
                zlabel(stmp.zname);
                colorbartitle(stmp.wname);
            else
                h = plot3c(stmp.x/xunits, stmp.y/yunits, stmp.z/zunits, stmp.z/zunits, '.');
                zlabel(stmp.zname);
                colorbartitle(stmp.zname);
            end
        end % plotz
            
        function zptilt = removetilt(S)
            
            % remove the holes
            iuse = ~isnan(S.z);
            [xp, yp, zp] = filterdata(...
                iuse,...
                S.x, S.y, S.z);

            % calculate the best fit plane
            n = length(zp(:)); 
            A = [xp(:) yp(:) ones(n,1)];
            abc = A \ zp(:);
            
            % be sure not to change the holeval's when subtracting the
            % plane
            zptilt = zeros(size(S.z));
            zptilt(iuse) = zp(:) - A*abc;

            % add back the mean value for the constant offset
            zptilt(iuse) = zptilt(iuse) + mean(zp(:));
            
            % restore the holes
            zptilt(~iuse) = NaN;
            
        end % remove tilt
        
        function bResult = is4D(S)
            bResult = (length(S.w) == length(S.z));
        end %
        
        function S = SetHolevalNaN(S, holeval)
            % change points = holeval to NaN so that they will appear as
            % holes when plotted
            if ~isnan(holeval),
                if S.is4D,
                    [S.z, S.w] = filterdata(S.z ~= holeval, S.z, S.w, NaN);
                else
                    S.z = filterdata(S.z ~= holeval, S.z, NaN);
                end
            else
                if S.is4D,
                    S.w(isnan(S.z)) = NaN;
                end
            end
            
        end % SetHolevalNaN
        
        function S = ApplyUnits(S, Units)
            S.x = Units.xunits * S.x;
            S.y = Units.yunits * S.y;
            S.z = Units.zunits * S.z;
            
        end
        
        function S = FilterData(S, iuse)
            % in-place filtering
            if is4D(S),
                [S.x, S.y, S.z, S.w] = filterdata(iuse, S.x, S.y, S.z, S.w);
            else
                [S.x, S.y, S.z] = filterdata(iuse, S.x, S.y, S.z);
            end
            
        end % filterdata
        
        function Sout = cat(Sin)
            % Sin = array of CDataSet4D
            % e.g. Sout = cat([S1 S2])
            Sout = CDataSet4D(...
                cat(1,Sin(:).x) ...
                , cat(1,Sin(:).y) ...
                , cat(1,Sin(:).z) ...
                , Sin(1).zname ...
                , cat(1,Sin(:).w) ...
                , Sin(1).wname ...
                );

        end % cat
        
        function S = OutlierFilter(S, z_or_w, n_sig)
        
            vtmp = S.(z_or_w);
            itmp = (1:length(vtmp))';
            if length(vtmp) <= 3, return, end
            [p, s, mu] = polyfit(itmp, vtmp, 2);
            vtmp = vtmp - polyval(p, itmp, [], mu);
            vmu = mean(vtmp);
            vsig = std(vtmp);
            S = S.FilterData(vtmp > vmu - n_sig*vsig & vtmp < vmu + n_sig*vsig);
                    
        end % OutlierFilter
        
        function S = Resample(S, b, N)
            % [b, a] = filter coefficients, using filtfilt()
            % filter S.z, then decimate each field
            % Sout.x = Sin.x(1:N:end)
            %
            
            %S = S.filtfilt(b,a);
            %S.z = resample(S.z, 1, N);
            S.z = DelayCompensatedFIRFilter(b, S.z);

            % decimate
            if is4D(S), S.w = S.w(1:N:end); end
            S.x = S.x(1:N:end);
            S.y = S.y(1:N:end);
            S.z = S.z(1:N:end);
                             
        end % Resample
                        
        function S = filtfilt(S, b, a)
            % [b, a] = filter coefficients, using filtfilt()
            % filter S.z, then decimate each field
            % Sout.x = Sin.x(1:N:end)
            %
            if any(isnan(S.z(:))),
                warning('removing NaN values before filtering');
                S = S.FilterData(~isnan(S.z));
            end
            
            % remove best fit line
            xtmp = (0:length(S.z)-1)';
            [p, s, mu] = polyfit(xtmp, S.z, 1);
            zline = polyval(p, xtmp, [], mu);
            ztmp = S.z - zline;
            
            % apply filter
            if length(ztmp) > 3 * (max(length(b),length(a))+1),
                ztmp = filtfilt(b, a, ztmp);
            else
                % do nothing
            end
            
            % add line back
            S.z = ztmp + zline;
                        
        end % Resample

        function Sdiff = minus(S1, S2)
            % triangulate the more coarse grid
            if length(S1.x(:)) < length(S2.x(:)),
                F = TriScatteredInterp(S1.x, S1.y, S1.z);
                z1onS2 = F(S2.x, S2.y);
                
                Sdiff = CDataSet4D(S2.x, S2.y, z1onS2 - S2.z, [S1.zname ' - ' S2.zname]);
                
            else
                F = TriScatteredInterp(S2.x, S2.y, S2.z);
                z2onS1 = F(S1.x, S1.y);
                
                Sdiff = CDataSet4D(S1.x, S1.y, S1.z - z2onS1, [S1.zname ' - ' S2.zname]);
                
            end
        end % minus
        
    end % end methods
    
end % classdef

    