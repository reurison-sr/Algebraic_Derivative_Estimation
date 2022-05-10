
%{
"A derivative estimation toolbox based on algebraic
methods theory and practice" - 2007. 
Authors: Josef Zehetner, Johann Reger and Martin Horn
%}

clear all; close all; clc

syms tau T
N = 2; % N : order off the Taylor-series expasion.
j = 1; % j : order of derivative for the input signal
v = 0; % v : number of additional integrals.

%Equation (11)
PIjNv = 0;
for k1=0:1:(N-j)
    for k2=0:1:j
        num = ((T - tau)^(v + k1 + k2))*(-tau)^(N - k1 - k2);
        den = factorial(k1)*factorial(k2)*factorial(N-j-k1)*factorial(j-k2)*...
            factorial(N-k1-k2)*factorial(v+k1+k2)*(N-k1+1);
        PIjNv = PIjNv + (num/den);
    end
end

q = ( factorial(N+j+v+1)*factorial(N+1)*(-1)^j ) / (T^(N+j+v+1));
aj = ((-1)^j)*q*Sum