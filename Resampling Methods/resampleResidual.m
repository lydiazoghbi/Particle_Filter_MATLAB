function [ indx ] = resampleResidual(w, N)
% [ indx ] = resampleResidual(w, N)
% Residual resampling method for particle filtering.  
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
% Output:
%       indx the resampled index according to the weight sequence

if nargin == 1
   N = length(w);
end
M = length(w);
w = w / sum(w);
indx = zeros(1, N);
% integer parts:
Ns = floor(N.* w);
R = sum(Ns);
% Draw the deterministic part:
i = 1;
j = 0;
while j < M
    j = j + 1;
  cnt = 1;
  while cnt <= Ns(j)
    indx(i) = j;
    i = i + 1; cnt = cnt + 1;
  end;
end;
% The fractions: Multinomial) resampling
N_rdn = N - R;
Ws =(N*w - Ns)/N_rdn;
Q = cumsum(Ws);
while(i <= N)
    sampl = rand;  %(0,1]
    j = 1;
    while(Q(j) < sampl),
        j = j + 1;
    end;
    indx(i) = j;
    i = i + 1;
end

