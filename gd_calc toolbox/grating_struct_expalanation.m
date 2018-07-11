grating.pmt{1:NumberOfPermittivities} % a look up table of permittivities
grating.pmt_sub_index % index into grating.pmt{} for substrate
grating.pmt_sup_index % index into grating.pmt{} for supertrate
% coefficients for directional grating vectors
[grating.d21 grating.d31] % vector d1
[grating.d22 grating.d32] % vector d2
grating.stratum{1:N} % numbered from the bottom up

% struct stratum:
% for homogeneous
stratum.type = 0; % homogeneous
stratum.thick
stratum.pmt_index % index into grating.pmt{}

% for uniperiodic
stratum.type = 1;
stratum.thick
[stratum.h11; stratum.h12] % = [ [d21; d31]./d1^2 [d22; d32]./d2^2 ] \ [d_stratum/dstratum^2]
stratum.stripe{1:Nstripes}
% a stripe for uniperiodic:
stripe.c1 % fraction of d_stratum
stripe.pmt_index % index into grating.pmt{}

