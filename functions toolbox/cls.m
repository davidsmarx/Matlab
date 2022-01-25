function cls(varargin)
% cls
% close all figure windows and clear all variables
% cls hidden closes hidden windows also
%
% CheckOption('except', {}, varargin{:}); % cell array of 'var' not to clear

list_except = CheckOption('except', {}, varargin{:});

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

if ~isempty(list_except),
    clearvars('-except', list_except{:});
else
    evalin('base','clear');
end

evalin('base','constants');

evalin('base','unitsdefinitions');