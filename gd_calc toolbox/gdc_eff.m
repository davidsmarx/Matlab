function [R,T]=gdc_eff(scat_field,inc_field)
% gdc_eff
%
% Compute grating diffraction efficiencies from gdc output. (See gdc.m
% comment header.)
%
% Syntax:
%
%   [R,T]=gdc_eff(scat_field,inc_field);
%
% Note: In the following output specification, s and p are the incident
% field's polarization basis vectors (GD-Calc.pdf, equations 4.16-19). For
% an incident E field of arbitrary complex amplitude A*s+B*p, the
% diffraction efficiency in a particular diffraction order (reflected or
% transmitted) is
%   (abs(A)^2*eff1 ...
%   +abs(B)^2*eff2 ...
%   +real(conj(A)*B)*(2*eff3-eff1-eff2) ...
%   -imag(conj(A)*B)*(2*eff4-eff1-eff2)) ...
%   /(abs(A)^2+abs(B)^2)
% (eff1 ... eff4 correspond to the data fields of R or T for a particular
% diffraction order.)
%
% inputs:
%
%   scat_field,inc_field: from gdc
%
% outputs:
%
%   R (size-[1,?] struct): data for reflection efficiencies
%
%     R(k).m1, m2 (integer): diffraction order indices
%
%     R(k).eff1 (parameter, real): efficiency for incident field amplitude
%     = s (i.e., A=1, B=0; linear polarization, TE)
%
%     R(k).eff2 (parameter, real): efficiency for incident field amplitude
%     = p (i.e., A=0, B=1; linear polarization, TM)
%
%     R(k).eff3 (parameter, real): efficiency for incident field amplitude
%     = (s+p)/sqrt(2) (i.e., A=B=1/sqrt(2); linear polarization, diagonal)
%
%     R(k).eff4 (parameter, real): efficiency for incident field amplitude
%     = (s-i*p)/sqrt(2) (i.e., A=1/sqrt(2), B=-i/sqrt(2); right-hand
%     circular polarization)
%
%   T (size-[1,?] struct): data for transmission efficiencies
%
%     T.m1, m2, eff1, eff2, eff3, eff4: same format as R
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

len=length(scat_field);
R=struct('m1',{scat_field.m1},'m2',{scat_field.m2});
T=R;
for k=1:len
    [f1i,f1r,f1t,Rss,Rsp,Rps,Rpp,Tss,Tsp,Tps,Tpp]=...
        repmat_extend(inc_field.f1,scat_field(k).f1r,scat_field(k).f1t,...
        scat_field(k).Rss,scat_field(k).Rsp,...
        scat_field(k).Rps,scat_field(k).Rpp,...
        scat_field(k).Tss,scat_field(k).Tsp,...
        scat_field(k).Tps,scat_field(k).Tpp);
    % GD-Calc.pdf, equations 4.38, 4.33
    R(k).eff1=abs(Rss).^2+abs(Rps).^2;
    R(k).eff2=abs(Rsp).^2+abs(Rpp).^2;
    R(k).eff3=0.5*(abs(Rss+Rsp).^2+abs(Rps+Rpp).^2);
    R(k).eff4=0.5*(abs(Rss-i*Rsp).^2+abs(Rps-i*Rpp).^2);
    c=-real(f1r)./real(f1i);
    R(k).eff1=R(k).eff1.*c;
    R(k).eff2=R(k).eff2.*c;
    R(k).eff3=R(k).eff3.*c;
    R(k).eff4=R(k).eff4.*c;
    % GD-Calc.pdf, equations 4.39, 4.34
    T(k).eff1=abs(Tss).^2+abs(Tps).^2;
    T(k).eff2=abs(Tsp).^2+abs(Tpp).^2;
    T(k).eff3=0.5*(abs(Tss+Tsp).^2+abs(Tps+Tpp).^2);
    T(k).eff4=0.5*(abs(Tss-i*Tsp).^2+abs(Tps-i*Tpp).^2);
    c=real(f1t)./real(f1i);
    T(k).eff1=T(k).eff1.*c;
    T(k).eff2=T(k).eff2.*c;
    T(k).eff3=T(k).eff3.*c;
    T(k).eff4=T(k).eff4.*c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [varargout]=repmat_extend(varargin)
% [A1,A2,...]=repmat_extend(A1,A2,...);
% Repmat-extend arguments' singleton dimensions to match array sizes.
% (Scalars will not be repmat-extended.) Size compatibility of input
% arguments is assumed.
n=nargin; % >=nargout
s=cell(1,n);
s_=[1,1];
for j=1:n
    s{j}=size(varargin{j});
    len=length(s{j});
    s_(end+1:len)=1;
    s_(1:len)=max(s_(1:len),s{j});
end
ndims_=length(s_);
n=nargout;
varargout=cell(1,n);
for j=1:n
    if all(s{j}==1)
        varargout{j}=varargin{j}; % (Do not repmat-extend scalars.)
    else
        s{j}(end+1:ndims_)=1;
        ext=s_./s{j};
        if all(ext==1)
            varargout{j}=varargin{j};
        else
            varargout{j}=repmat(varargin{j},ext);
        end
    end
end
