function [ X ] = spectralEstCVX( meas, cameraMat, basisFcns, alpha, beta, gamma )

nWaves = size(basisFcns,1);
nBasis = size(basisFcns,2);

h = size(meas,1);
w = size(meas,2);
nChannels = size(meas,3);

meas = reshape(meas,h*w,nChannels)';



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


b = [meas; zeros(nWaves-1,h*w)];
A = [cameraMat; sqrt(alpha)*Rlambda]*basisFcns;

cvx_begin
    variable X(nBasis,h*w)
    if beta ~= 0
        minimize sum(sum_square(b - A*X)) + beta*sum(sum(abs(RX*(X)'))) + beta*sum(sum(abs(RY*(X)')))
    else
        minimize sum(sum_square(b - A*X))    
    end
    subject to
        if gamma > 0
            0 <= basisFcns*X <= 1
        end 
cvx_end


end

