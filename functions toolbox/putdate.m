% PUTDATE and a mouse click on a figure adds a time and date
% stamp in the desired location using gtext.
% Location is very predicatble on 2D plots,
% less so on 3D plots.
% Like datestamp.m by Peter Newbury except that
% 1) it has been tested in MATLAB 5.1.0.421 and MATLAB 4.
% 2) it puts the date stamp anywhere on the plot.

% Mariano Garcia 12/5/96 garcia@tam.cornell.edu
% please send bug info to: garcia@tam.cornell.edu

% get time and date info string from MATLAB and
% convert it to a more appropriate text string

datetime=fix(clock);
month=datetime(2);
if month ==1 monthstr='Jan'; end;
if month==2 monthstr='Feb'; end;
if month==3 monthstr='Mar'; end;
if month==4 monthstr='Apr'; end;
if month==5 monthstr='May'; end;
if month==6 monthstr='Jun'; end;
if month==7 monthstr='Jul'; end;
if month==8 monthstr='Aug'; end;
if month==9 monthstr='Sep'; end;
if month==10 monthstr='Oct'; end;
if month==11 monthstr='Nov'; end;
if month==12 monthstr='Dec'; end;

ampmstr='AM';
yearstr=int2str(datetime(1));
daystr=num2str(datetime(3)/100);
% take care of 10, 20, and 30. First two elements in string are "0."
if length(daystr)<4
  daystr(4)='0'; end;
daystr=daystr(3:4);

% adjust AM-PM and 12AM.
hourstr=int2str(datetime(4));
if datetime(4) > 12
   hourstr = int2str(datetime(4)-12);
   ampmstr='PM';
end;
if datetime(4)==0
  hourstr = '12';
end;

% take care of minutes. "00" is a possibility.
minstr=num2str(datetime(5)/100);
if length(minstr)<4
  minstr(4)='0';
   if datetime(5)==0 minstr='0000'; end;
end;

minstr=minstr(3:4);

colon=':';
space=' ';

[x,y]=ginput(1);
datestr=str2mat(daystr,monthstr,yearstr,space,...
hourstr,colon,minstr,ampmstr);
datestr=datestr';
datestr=datestr(:)';

% get rid of some spaces
for i=22:29
datestr(i)=datestr(i+3);
end;
for i=19:30
datestr(i)=datestr(i+2);
end;
for i=13:30
datestr(i)=datestr(i+2);
end;
datestr=datestr(1:25);

% print the thing on the screen
output=str2mat(datestr);
text(x,y,output);
