function [ indx ] = resampleMetropolisPP( w, B )
% [ indx ] = resampleMetropolisPP( w, B )
% Metropolis (parallel) resampling method for particle filtering, 
% This method is designed for parallel processing but here IT runs in serial.  
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       B    the number of iterations to be performed for each particle before settling on its chosen ancestor
% Output:
%       indx  the resampled index according to the weight sequence

if nargin == 1
  B = length(w) / 10;
end
M = length(w);

U = 0.5 + (M-0.5).*rand(M*B,1);
Indj = round(U);
indx = zeros(1, M);

b = 1;
j = 0;
while j < M
    j = j + 1;
    k = j;
    while b <= B
        i = Indj(b);
        b = b + 1;
        if rand <= w(i)/w(k)
            k = i;
        end
    end
    indx(j) = k;
end
