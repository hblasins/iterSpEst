function [ X ] = spectralEstCVX( meas, cameraMat, basisFcns, alpha, beta, gamma )

% [ X ] = spectralEstCVX( meas, cameraMat, basisFcns, alpha, beta, gamma )
%
% This function computes the reflectance basis function weights that best
% explain the data in meas using cvx toolbox (). This function solves the same
% optimization problem as spectralEstADMM. However it will work only for
% small problems, where the total number of variables is about 50,000. This
% function was created for comparison purposes only. More details about the
% algorithm can be found in
% 
% H. Blasinski, J. Farrell and B. Wandell; 'Iterative Algorithm for
% spectral estimation with spatial smoothing,' ICIP2015, Quebec City.
%
%
% Parameters:
%  meas - a h x w x c 3D matrix, where the first two dimensions represent
%      image height and width. The last dimension is the number of imaging
%      channels
%
%  cameraMat - a c x nWaves matrix where each row represents the spectral
%      responsivity of a given camera channel. It is assumed that the
%      spectrum is discretized to nWaves bins.
%
%  basisFcns - a nWaves x nBasis matrix where each column is a reflectance
%      basis function discretized to nWaves spectral bins.
%
%  alpha - a scalar controling the spectral smoothess of the solution.
%  
%  beta - a scalar controlling the spatial smoothness of the solution.
%
%  gamma - a scalar controlling whether the solution is bounded to the
%     [0,1] interval. If gamma is equal to zero this constraint is NOT
%     enforced. If gamma is different from zero then the constraint IS
%     enfoced and the value of the parameter does not matter.
%  
%
% Returns:
%  X - a h x w x nBasis matrix of estimated spectral reflectance weights.
%  
%
% Copyright, Henryk Blasinski 2015.


nWaves = size(basisFcns,1);
nBasis = size(basisFcns,2);

h = size(meas,1);
w = size(meas,2);
nChannels = size(meas,3);

meas = reshape(meas,h*w,nChannels)';


% Spectral gradient
Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];


% y direction spatial gradient
Ry = [eye(h-1) zeros(h-1,1)] - [zeros(h-1,1) eye(h-1)];
RY = sparse([]);
for i=1:w;
    RY = blkdiag(RY,Ry);
end

% x direction spatial gradient
Rx = zeros(h,2*h);
for y=1:h
   Rx(y,:) = [zeros(1,y-1) 1 zeros(1,h-y) zeros(1,y-1) -1 zeros(1,h-y)];
end

RX = sparse([]);
for x=1:(w-1)
   RX((x-1)*h+1:x*h,(x-1)*h+1:(x+1)*h) = Rx; 
end


cvx_begin
    variable X(nBasis,h*w)
    if beta ~= 0
        minimize sum(sum_square(meas - cameraMat*basisFcns*X)) + alpha*sum(sum_square(Rlambda*basisFcns*X)) + beta*sum(sum(abs(RX*(X)'))) + beta*sum(sum(abs(RY*(X)')))
    else
        minimize sum(sum_square(meas - cameraMat*basisFcns*X)) + alpha*sum(sum_square(Rlambda*basisFcns*X))   
    end
    subject to
        if gamma > 0
            0 <= basisFcns*X <= 1
        end 
cvx_end

X = reshape(X',[h w nBasis]);

end

