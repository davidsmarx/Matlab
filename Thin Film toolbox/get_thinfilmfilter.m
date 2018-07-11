function [d, n] = get_thinfilmfilter(ftype,refl)
% [d, n] = get_thinfilmfilter(filtertype,refl)
% valid filtertype:
%     'Natick100GHz' (0)
%     'Marlborough200GHz' (1)
%     'Keqi50GHz' (2)
%     'Yadlowsky50GHz'
% refl = reference wavelength in vacuum
%
% d = layer thicknesses [same units as refl, d = 0.25*refl/real(n)]
% n = layer indexes
% order: substrate medium is first, incident medium is last
% length(n) = length(d) + 2;

if nargin < 2,
   refl = 1.55; % um reference wavelength
   if nargin == 0,
      ftype = 0;
   end
end
fprintf('get filter: ref lam = %f\n',refl);

switch ftype,
   
case {0, 'Natick100GHz'}, % Natick 100 GHz
   % sub /  
   % (HL)^7 H 8L H (LH)^7 L  
   % [(HL)^8 H 8L H (LH)^8 L] ^2  
   % (HL)^7 H 8L H (LH)^6 L  0.56H 1.2 L  
   % / air
   % where nsub = 1.66  nH = 2.07   nL = 1.465 0.7106H 0.5618L
   nh = 2.07; nl = 1.465;
   d0 = [1 1]';
   d1 = [repmat(d0,7,1); 1; 8; 1; repmat(d0,7,1); 1];
   n1 = repmat([nh; nl],16,1);
   d2 = [repmat(d0,8,1); 1; 8; 1; repmat(d0,8,1); 1];
   n2 = repmat([nh; nl],18,1);
   d3 = [repmat(d0,8,1); 1; 8; 1; repmat(d0,8,1); 1];
   n3 = repmat([nh; nl],18,1);
   d4 = [repmat(d0,7,1); 1; 8; 1; repmat(d0,6,1); 1; 0.6; 0.7];
   n4 = repmat([nh; nl],16,1);
   
   nt = 1.66;  % index of transmission medium

case {1, 'Marlborough200GHz'}, % Marlborough 200Ghz
   %( (1H 1L)4  1H 16L 1H (1L 1H)4 1L)
   %( (1H 1L)5  1H 12L 1H (1L 1H)5 1L)
   %( (1H 1L)5  1H 12L 1H (1L 1H)5 1L) 
   %( (1H 1L)4  1H 16L 1H (1L 1H)4 1L)
   % 
   % 1.717H 2.683L
   % 
   %where H means one quarter wave optical thickness of Nb2O5, n=2.26
   %and     L means one quarter wave optical thickness of SiO2, n = 1.45
   nh = 2.26; nl = 1.45;
   d0 = [1 1]';
   d1 = [repmat(d0,4,1); 1; 16; 1; repmat(d0,4,1); 1];
   n1 = repmat([nh; nl],10,1);
   d2 = [repmat(d0,5,1); 1; 12; 1; repmat(d0,5,1); 1];
   n2 = repmat([nh; nl],12,1);
   d3 = [repmat(d0,5,1); 1; 12; 1; repmat(d0,5,1); 1];
   n3 = repmat([nh; nl],12,1);
   d4 = [repmat(d0,4,1); 1; 16; 1; repmat(d0,4,1); 1; 1.717; 2.683];
   n4 = repmat([nh; nl],11,1);
   
   nt = 1.66;  % index of transmission medium

case {2, 'Keqi50GHz'},
   % (1H 1L)8 1H 8L 1H (1L 1H)8 1L 
   % (1H 1L)9 1H 8L 1H (1L 1H)9 1L 
   % (1H 1L)9 1H 8L 1H (1L 1H)9 1L 
   % (1H 1L)8 1H 8L 1H (1L 1H)7 1L .6H .7L
   % H = 2.07
   % L = 1.465
   % Sub = 1.66
   nh = 2.07; nl = 1.465;
   d0 = [1 1]';
   d1 = [repmat(d0,8,1); 1; 8; 1; repmat(d0,8,1); 1];
   n1 = repmat([nh; nl],length(d1)/2,1);
   d2 = [repmat(d0,9,1); 1; 8; 1; repmat(d0,9,1); 1];
   n2 = repmat([nh; nl],length(d2)/2,1);
   d3 = [repmat(d0,9,1); 1; 8; 1; repmat(d0,9,1); 1];
   n3 = repmat([nh; nl],length(d3)/2,1);
   d4 = [repmat(d0,8,1); 1; 8; 1; repmat(d0,7,1); 1; 0.6; 0.7];
   n4 = repmat([nh; nl],length(d4)/2,1);
   
   nt = 1.66;
   
case {3, 'Yadlowsky50GHz'},
   nh = 2.0655; nl = 1.457;
   
   d0 = [1 1]';
   d1 = [repmat(d0,8,1); 1; 8; 1; repmat(d0,8,1); 1];
   n1 = repmat([nh; nl],length(d1)/2,1);
   d2 = [repmat(d0,9,1); 1; 6; 1; repmat(d0,9,1); 1];
   n2 = repmat([nh; nl],length(d2)/2,1);
   d3 = [repmat(d0,9,1); 1; 6; 1; repmat(d0,8,1); 1];
   n3 = repmat([nh; nl],length(d3)/2,1);
   d4 = [repmat(d0,9,1); 1; 8; 1; repmat(d0,7,1); 1; 1.37415; 1.28443];
   n4 = repmat([nh; nl],length(d4)/2,1);
   
   nt = 1.66;
otherwise,
   error('unknown filter type');
end

d = [d1; d2; d3; d4]; 
nr = [n1; n2; n3; n4]; ni = 0;

% convert thicknesses from units of quarter-waves into um
d = 0.25*refl*(d./nr);

n = [nt; nr + j*ni; 1];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other filter definitions
% the first definition is identical to the one above.
if 0,  % get filter prescription from file
   [str, N] = textread('dfc_filter.txt','%s %d',1,'headerlines',4);
   fprintf('filter has %d layers\n',N);
   
   [lth, d, lin, nr, ni] = textread('dfc_filter.txt','%s %f %s %f %f','headerlines',5);
   
   % reverse order of layers from substrate-to-incident medium to incident medium-to-substrate
   d = d(end:-1:1); nr = nr(end:-1:1); ni = ni(end:-1:1);
   
else,  % use the following prescription
   
   % ( (1H 1L)^5  3H 12L 1H (1L 1H)^5 1L)^1
   % ( (1H 1L)^6  1H 6L 1H (1L 1H)^6 1L)^1
   % ( (1H 1L)^6  1H 10L 1H 1L 3H (1L 1H)^5 1L)^1 
   % ( (1H 1L)^5  1H 10L 1H (1L 1H)^6 1L  )^1
   % where H is one quarter wave of high index n=2.26
   % L is one quarterwave of low index n= 1.46
   nh = 2.26; nl = 1.46;
   d0 = [1 1]';
   d1 = [repmat(d0,5,1); 3; 12; 1; repmat(d0,5,1); 1];
   n1 = repmat([nh; nl],12,1);
   d2 = [repmat(d0,6,1); 1; 6; 1; repmat(d0,6,1); 1];
   n2 = repmat([nh; nl],14,1);
   d3 = [repmat(d0,6,1); 1; 10; 1; 1; 3; repmat(d0,5,1); 1];
   n3 = repmat([nh; nl],14,1);
   d4 = [repmat(d0,5,1); 1; 10; 1; repmat(d0,6,1); 1];
   n4 = repmat([nh; nl],13,1);
   
end

