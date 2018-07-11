function cls(varargin)
% cls
% close all figure windows and clear all variables
% cls hidden closes hidden windows also

sss = whos('global');
if any(strcmp({sss.class},'activex')), 
   disp('Warning: There are Activex components open:');
   na = strmatch('activex',{sss.class});
   disp({sss(na).name});
   a = input('1) execute Exlclose(Excel);\n2) Quit cls\n3) continue cls\n');
   switch a,
   case 1, Exlclose(Excel);
   case 2, return,
   case 3,
   end
end

if strmatch('hidden',varargin),
    evalin('base','close all hidden');
else,
    evalin('base','close all');
end

evalin('base','clear');

evalin('base','constants');

evalin('base','unitsdefinitions');