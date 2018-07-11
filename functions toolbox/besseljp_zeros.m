function znk = besseljp_zeros(nu, k)

% http://mathworld.wolfram.com/BesselFunctionZeros.html
%
% The first few roots j_(n,k)^' of the derivative of the Bessel function
% J_n^'(x) are given in the following table for small nonnegative integer
% values of n and k. Versions of the Wolfram Language prior to 6 
% implemented these zeros as BesselJPrimeZeros[n, k] in the BesselZeros
% package which is now available for separate download (Wolfram Research).
% Note that contrary to Abramowitz and Stegun (1972, p. 370), the Wolfram
% Language defines the first zero of J_0^'(z) to be approximately 3.8317
% rather than zero. 
% 
% k J_0^'(x) J_1^'(x) J_2^'(x) J_3^'(x) J_4^'(x) J_5^'(x) 
% 1 3.8317 1.8412 3.0542 4.2012 5.3175 6.4156 
% 2 7.0156 5.3314 6.7061 8.0152 9.2824 10.5199 
% 3 10.1735 8.5363 9.9695 11.3459 12.6819 13.9872 
% 4 13.3237 11.7060 13.1704 14.5858 15.9641 17.3128 
% 5 16.4706 14.8636 16.3475 17.7887 19.1960 20.5755 

pkn = [
    3.8317 1.8412 3.0542 4.2012 5.3175 6.4156 
    7.0156 5.3314 6.7061 8.0152 9.2824 10.5199
    10.1735 8.5363 9.9695 11.3459 12.6819 13.9872
    13.3237 11.7060 13.1704 14.5858 15.9641 17.3128
    16.4706 14.8636 16.3475 17.7887 19.1960 20.5755
    ];

%besseljp_zeros = @(nu, k)(pkn(k, nu+1));

znk = pkn(k, nu+1);