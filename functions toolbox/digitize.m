function iXY = digitize()
%function iXY = digitize()
%                                                                            
%jdc 17-Feb-00
%PURPOSE   
%     > Interactively digitizes data from an image file.  
%     > X-Y axis scales can be any combination of linear, log.
%     > Automatic compensation for angular misalignment of  
%       image with respect to screen axes.
%INPUT
%       (All input is interactively prompted)
%     > imagename - image filespec (any Matlab-compatible format)
%     > three points defining coordinate axes  
%     > arbitrary number of points on graph - digitized with left mouse button  
%OUTPUT
%     > numbered sets of iXY triplets written to screen and optionally to disc storage
%     > iXY - 3 x i matrix of all points digitized from image 
%
%TYPICAL CALL 
%       digitize <RETURN>

if ~nargin,
   imagename = input('Enter the image filespec => ','s');
end   
pic = imread(imagename);
image(pic)
set(gcf,'Units','normalized','Position',[0 .15 1 .85])
set(gca,'Units','normalized','Position',[0   0 1   1]);
fprintf([ 'Select the plot origin with the left mouse button ...\n' 7])
[Xopixels,Yopixels] = ginput(1);
line(Xopixels,Yopixels,'Marker','.','Color','r')
OriginXYdata = input([   '  Enter the [x y] graph values at the origin => ' 7]);
if length(OriginXYdata)==1,
   OriginXYdata(2) = input(['       Enter the Y graph value at the origin => ' 7]);
end
fprintf(['\nSelect a point on the X-axis ...\n' 7])
[XAxisXpixels,XAxisYpixels] = ginput(1);
line(XAxisXpixels,XAxisYpixels,'Marker','.','Color','r')
XAxisXdata = input(['Enter the X-value of the point => ' 7]);
logx =     input(   ' X-axis scale? [ {lin} | log ] => ','s');
if findstr('log',lower(logx)),
   logx = 1;
   scalefactorXdata = log10(XAxisXdata/OriginXYdata(1));
else
   logx = 0;
   scalefactorXdata = XAxisXdata - OriginXYdata(1);
end 
th = atan((XAxisYpixels-Yopixels)/(XAxisXpixels-Xopixels));   % note image file line 1 is at top
%fprintf(          '              Graph rotation angle: %+3.1f degrees +cw\n', th*180/pi)
rotmat = [cos(th) sin(th); -sin(th) cos(th)];    % axis rotation matrix

fprintf(['\nSelect a point on the Y-axis ...\n' 7])
[YAxisXpixels,YAxisYpixels] = ginput(1);
line(YAxisXpixels,YAxisYpixels,'Marker','.','Color','r')
YAxisYdata = input(['Enter the Y-value of the point => ' 7]);
logy =     input(   ' Y-axis scale? [ {lin} | log ] => ','s');
if findstr('log',lower(logy)),
   logy = 1;
   scalefactorYdata = log10(YAxisYdata/OriginXYdata(2));
else
   logy = 0;
   scalefactorYdata = YAxisYdata - OriginXYdata(2);
end 
delxyx = rotmat*[(XAxisXpixels-Xopixels);(XAxisYpixels-Yopixels)];
delxyy = rotmat*[(YAxisXpixels-Xopixels);(YAxisYpixels-Yopixels)];
delXcal = delxyx(1);
delYcal = delxyy(2);

numberformat = '%6.2f';
nXY = [];
ng = 0;
while 1,
   fprintf(['\n Select graph points with left mouse button ...\n',...
              ' To terminate capture - select a point to the left of or below the graph axes\n' 7])
   n = 0;
   
%---------------------------------------------------------------------------------------- data acquire loop
   while 1
      [x,y] = ginput(1);                       
      line(x,y,'Marker','+','Color','g')
      xy = rotmat*[(x-Xopixels);(y-Yopixels)];
      delXpoint = xy(1);
      delYpoint = xy(2);
      if delXpoint>=0 & delYpoint<=0,      
         if logx,
            x = OriginXYdata(1)*10^(delXpoint/delXcal*scalefactorXdata);
         else
            x = OriginXYdata(1) + delXpoint/delXcal*scalefactorXdata;
         end
         if logy, 
            y = OriginXYdata(2)*10^(delYpoint/delYcal*scalefactorYdata);
         else  
            y = OriginXYdata(2) + delYpoint/delYcal*scalefactorYdata;
         end
         n = n+1;
         xpt(n) = x;
         ypt(n) = y;
         fprintf(      ['%4.0f  ' numberformat '  ' numberformat '\n'],n,x,y)
         ng = ng+1;
         nXY(ng,:) = [n x y];
      else
         break
      end
   end
   
%---------------------------------------------------------------------------------------------
   fname =        input([   '  Enter a filespec to save the data to (default: no save) => ' 7],'s');
   if ~isempty(fname),
      numberformat = input(['        Enter a format to save the data (default: %6.2f ) => ' 7],'s');
      if isempty(numberformat),
         numberformat = '%6.2f';
      end
      
      for i = 1:n,
         fprintf(fname,['%4.0f  'numberformat '  ' numberformat '\n'],i,xpt(i),ypt(i))
      end
      
      I = findstr('\',fname);
      if ~isempty(I),
         for i = length(I):-1:1,
            fname = [fname(1:i) '\' fname(i+1:length(fname))];
         end
      end   
      fprintf([   sprintf(['                         %3.0f x-y data pairs saved to file => ' fname],n) '\n'])
   end
   yn = input([            '          Digitize more data from this graph? [ {y} | n ] => ' 7],'s');
   if findstr('n',lower(yn)),
      break
   end
end   
