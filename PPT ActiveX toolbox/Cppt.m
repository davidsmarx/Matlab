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
        
        function [hPic, new_slide] = CopyFigNewSlide(S, hfig)
            % Add a new slide (with title object):
            slide_count = get(S.Presentation.Slides,'Count');
            slide_count = int32(double(slide_count)+1);

            layout = S.AppPpt.ActivePresentation.SlideMaster.CustomLayouts.Item(5); % was 10
            new_slide=S.AppPpt.ActivePresentation.Slides.AddSlide(slide_count, layout);

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
            
            options.Format = 'jpg';
            hgexport(hfig,'-clipboard', options);
            hPic = invoke(new_slide.Shapes,'Paste');

            %             % Set picture to width of slide, justify left, and 90% from the top, then close figure
            %             set(pic,'Width',720);
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

            %             % Set picture to width of slide, justify left, and 90% from the top, then close figure
            %             set(pic,'Width',720);
            %             set(pic,'Left',(slide_W-get(pic,'Width'))/2)
            %             set(pic,'Top',slide_H-slide_H*.905);
            
        end % CopyFigSlide

    end % methods
    
    
end % classdef