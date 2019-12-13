function dist = resampleMethod(w,N,method)
% This function is used to complete the resampling step of the particle
% filter method. w is the input weight vector for each of the particles.  N
% is the number of particles to generate in the new,unweighted
% distribution. "method" describes which resampling method will be used to
% generate the new sammples. 

if method==1
    dist = resampleMultinomial(w,N);
elseif method==2
    dist =  resampleResidual(w,N);
elseif method==3
    dist = resampleStratified(w,N);
elseif method==4
    dist =  resampleSystematic(w,N);
end