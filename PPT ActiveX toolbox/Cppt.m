classdef Cppt
    
    properties
        AppPpt
        Presentations
        Presentation
        listSlides
        
    end % properties
    
    
    methods
        function S = Cppt(fn, varargin)
            S.AppPpt = actxserver('PowerPoint.Application');
            S.AppPpt.Visible = 1;
            S.Presentations = S.AppPpt.Presentations;
            % invoke(S.Presentations)
            S.Presentation = S.Presentations.Add;
            % or open(S.Presentations, fn)
            
        end % Cppt
        
        function new_slide = NewSlide(S, slide_num)
            % new_slide = NewSlide(S)
            % new_slide = NewSlide(S, slide_num)
            %
            % if slide_num = [] or no input, add new slide at end

            layout = S.AppPpt.ActivePresentation.SlideMaster.CustomLayouts.Item(5); % was 10
            if ~exist('slide_num', 'var') || isempty(slide_num)
                slide_count = get(S.Presentation.Slides,'Count');
                slide_num = int32(double(slide_count)+1);
            end % if isempty(slide_num)

            new_slide=S.AppPpt.ActivePresentation.Slides.AddSlide(slide_num, layout);

        end % NewSlide
        
        function [hPic, new_slide] = CopyFigNewSlide(S, hfig)
            % [hPic, new_slide] = CopyFigNewSlide(S, hfig)
            %
            % add a new slide at end and copy fig to it

            % Add a new slide (with title object):
            new_slide = S.NewSlide([]);

            % Get height and width of slide:
            slide_H = S.Presentation.PageSetup.SlideHeight;
            slide_W = S.Presentation.PageSetup.SlideWidth;
            
            % Insert text into the title object:
            %             titletext = get(figure(nfigs(n)),'name');
            %             set(new_slide.Shapes.Title.TextFrame.TextRange,'Text',titletext);
            %             set(new_slide.HeadersFooters.Footer,'Visible',true)
            %             set(new_slide.HeadersFooters.Footer,'Text','|  Document Name  |  Enter Date')
            %             set(new_slide.HeadersFooters.SlideNumber,'Visible',true)
            %             set(nfigs(n),'units','normalized','outerposition',[0 0 1 1]);

            % from the Excel method:
            % figure(hfig); drawnow;
            % set(hfig,'invertHardcopy','off')
            % print -dbitmap
            %
            % rh = sheet.Range(cellid);
            % sheet.Paste(rh);

            % Copy figure(s) as bitmap(-dbitmap) (use '-dmeta' to copy as meta file instead)
            figure(hfig); drawnow;
            set(hfig,'invertHardcopy','off')
            
            %options.Format = 'jpg';
            copygraphics(hfig, 'ContentType', 'image'); %,'-clipboard', options);
            hPic = invoke(new_slide.Shapes,'Paste');

            %             % Set picture to width of slide, justify left, and 90% from the top, then close figure
            %             set(pic,'Width',720);
            set(hPic, 'Height', slide_H); % makes height of pic = height of slide
            if hPic.Width > slide_W,
                % then width is the limit, not height
                set(hPic, 'Width', slide_W)
                set(hPic, 'Left', 0)
            else
                % if height is the limit:
                set(hPic, 'Top', 0); % puts top of pic at top of slide
            end
            
            %             set(pic,'Left',(slide_W-get(pic,'Width'))/2)
            %             set(pic,'Top',slide_H-slide_H*.905);
            
        end % CopyFigNewSlide

        function hPic = CopyFigSlide(S, slide, hfig)
            % slide can be new_slide returned by CopyFigNewSlide()
            % slide can be S.AppPpt.ActivePresentation.Slides.Item(slidenum)
            
            % check that slide number slidenum exists
            %             slide_count = get(S.Presentation.Slides,'Count');
            %             if slidenum > slide_count
            %                 error(['only ' num2str(slide_count) ' slides in presentation']);
            %             end
            %
            %             slide = S.AppPpt.ActivePresentation.Slides.Item(slidenum);
            
            %             % Get height and width of slide:
            %             slide_H = S.Presentation.PageSetup.SlideHeight;
            %             slide_W = S.Presentation.PageSetup.SlideWidth;

            % Insert text into the title object:
            %             titletext = get(figure(nfigs(n)),'name');
            %             set(new_slide.Shapes.Title.TextFrame.TextRange,'Text',titletext);
            %             set(new_slide.HeadersFooters.Footer,'Visible',true)
            %             set(new_slide.HeadersFooters.Footer,'Text','|  Document Name  |  Enter Date')
            %             set(new_slide.HeadersFooters.SlideNumber,'Visible',true)
            %             set(nfigs(n),'units','normalized','outerposition',[0 0 1 1]);

            % Copy figure(s) as bitmap(-dbitmap) (use '-dmeta' to copy as meta file instead)
            figure(hfig); drawnow;
            set(hfig,'invertHardcopy','off')
            
            options.Format = 'jpg';
            hgexport(hfig,'-clipboard', options);
            hPic = invoke(slide.Shapes,'Paste');

            % Get height and width of slide:
            slide_H = S.Presentation.PageSetup.SlideHeight;
            slide_W = S.Presentation.PageSetup.SlideWidth;
            %             % Set picture to width of slide, justify left, and 90% from the top, then close figure
            %             set(pic,'Width',720);
            %             set(pic,'Left',(slide_W-get(pic,'Width'))/2)
            %             set(pic,'Top',slide_H-slide_H*.905);

            set(hPic, 'Height', slide_H); % makes height of pic = height of slide
            set(hPic, 'Top', 0); % puts top of pic at top of slide
            hPic.Left = 0;

        end % CopyFigSlide

        function hPic = AddPictureNewSlide(S, fn)
            % hPic = AddPictureNewSlide(S, fn)
            %
            % fn must be full path

            slide = S.NewSlide([]); % add new slide at end
            % Get height and width of slide:
            slide_H = S.Presentation.PageSetup.SlideHeight;
            slide_W = S.Presentation.PageSetup.SlideWidth;

            hPic = invoke(slide.Shapes, 'AddPicture', fn, true, true, 100, 100);
            % bring to front
            hPic.ZOrder(0);

            % scale the pic
            scale_w = slide_W / hPic.Width;
            scale_h = slide_H / hPic.Height;

            if scale_w < scale_h
                % width is the limit
                % apparently, scaling width also scales height
                hPic.Width = scale_w * hPic.Width;
                hPic.Left = 0;
                % center vertical
                hPic.Top = 0.5*(slide_H - hPic.Height);
            else
                % height is the limit
                % apparently, scaling width also scales height
                hPic.Width = scale_h * hPic.Width;
                hPic.Top = 0; % puts top of pic at top of slide
                % center horizontal
                hPic.Left = 0.5*(slide_W - hPic.Width);
            end

        end % AddPicture

    end % methods
    
    
end % classdef