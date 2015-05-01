function [ X, hist ] = spectralEstADMM( meas, cameraMat, basisFcns, alpha, beta, gamma, varargin )


p = inputParser;
p.addParamValue('maxIter',100,@isscalar);
p.addParamValue('tol',1e-3,@isscalar);
p.addParamValue('mu',10,@isscalar);
p.addParamValue('tauIncr',5,@isscalar);
p.addParamValue('tauDecr',5,@isscalar);
p.addParamValue('rhoInit',0.1,@isscalar);
p.addParamValue('rescaleRho',true,@islogical);
p.addParamValue('reference',[]);
p.addParamValue('tolConv',0,@isscalar);
p.addParamValue('verbose',false);
p.parse(varargin{:});
inputs = p.Results;

nBasis = size(basisFcns,2);
nWaves = size(basisFcns,1);
h = size(meas,1);
w = size(meas,2);
nChannels = size(meas,3);

meas = reshape(meas,h*w,nChannels)';

Z1 = zeros(nBasis,(h-1)*w + (w-1)*h);
Z1minus = zeros(size(Z1));
U1 = zeros(nBasis,(h-1)*w + (w-1)*h);

Z2 = zeros(nWaves,h*w);
Z2minus = zeros(size(Z2));
U2 = zeros(nWaves,h*w);

hist.prRes = zeros(inputs.maxIter,1);
hist.dualRes = zeros(inputs.maxIter,1);
hist.rho = inputs.rhoInit*ones(inputs.maxIter+1,1);
hist.pcg.iter = zeros(inputs.maxIter,1);
hist.pcg.acc = zeros(inputs.maxIter,1);


if isempty(inputs.reference) == false
    hist.conv = zeros(inputs.maxIter,1);
end

X = [];

for i=1:inputs.maxIter
    t1 = tic;
      
    
    % X update using cg
    t2 = tic;
    Atb = applyAt(meas,Z1-U1,Z2-U2,h,w,cameraMat,basisFcns, beta, gamma, hist.rho(i)/2);
    AtAhndl = @(x) applyAtA(x, h, w, cameraMat, basisFcns, alpha, beta, gamma, hist.rho(i)/2);
    [X, ~, hist.pcg.acc(i), hist.pcg.iter(i)] = pcg(AtAhndl,Atb,1e-6,10000,[],[],X(:));
    t3 = toc(t2);
    
    X = reshape(X,nBasis,h*w);
    msImg = reshape(X',[h,w,nBasis]);

    
    
    % Z1 update (spatial smoothness)
    
    tmpX = imfilter(msImg,[1 -1 0]);
    tmpY = imfilter(msImg,[1 -1 0]');
    
    
    dX = [reshape(tmpX(:,2:w,:),h*(w-1),nBasis)' reshape(tmpY(2:h,:,:),(h-1)*w,nBasis)'];
    
    % Anisotropic TV: Soft thresholding on Z1
    tmp = U1 + dX;
    if beta ~= 0
        Z1 = sign(tmp).*max(abs(tmp) - beta/hist.rho(i),0);
    else
        Z1 = tmp;
    end
    
    % Box constraint: Projection on Z2
    bFX = basisFcns*X;
    tmp = bFX + U2;
    if gamma ~= 0
        Z2 = max(min(tmp,1),0);
    else
        Z2 = tmp;
    end
    
    
    % Scaled dual variable update
    res1 = dX - Z1;
    res2 = bFX - Z2;
    U1 = U1 + res1;
    U2 = U2 + res2;
    
    
    % Residual computation
    hist.prRes(i) = sqrt(norm(res1,'fro')^2 + norm(res2,'fro')^2);
    hist.dualRes(i) = hist.rho(i)*sqrt(norm(Z1-Z1minus,'fro')^2 + norm(Z2-Z2minus,'fro')^2);
    
    
    % Re-scale the parameter rho
    if hist.prRes(i) > inputs.mu*hist.dualRes(i) && inputs.rescaleRho == true
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    end;
    if hist.dualRes(i) > inputs.mu*hist.prRes(i) && inputs.rescaleRho == true
        hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
    end;
    
    if inputs.verbose == true
        fprintf('ADMM iter %i (%f), primal res %e, dual res %e \n',i,toc(t1),hist.prRes(i),hist.dualRes(i));
        fprintf('     -> PCG (%f), err %f, nIter %i\n',t3,hist.pcg.acc(i),hist.pcg.iter(i));
    end
    
    if max(hist.prRes(i),hist.dualRes(i)) < inputs.tol
        break;
    end;
    
    if isempty(inputs.reference) == false
        hist.conv(i) = norm(X(:) - inputs.reference(:));
        if inputs.verbose == true
            fprintf('     -> True error %f (%f)\n',hist.conv(i),inputs.tolConv);
        end
        if hist.conv(i) < inputs.tolConv, break; end;
    end
    
    % Now we need to re-scale the scaled dual variable U as well as the
    % dual residuals
    U1 = U1*hist.rho(i)/hist.rho(i+1);
    U2 = U2*hist.rho(i)/hist.rho(i+1);
    
    Z1minus = Z1;
    Z2minus = Z2;
    
end

% If we terminate earlier, remove the unused portions of the vectors.
hist.prRes = hist.prRes(1:i);
hist.dualRes = hist.dualRes(1:i);
hist.rho = hist.rho(1:i);
hist.pcg.Iter = hist.pcg.iter(1:i);
hist.pcg.acc = hist.pcg.acc(1:i);

if isempty(inputs.reference) == false
    hist.conv = hist.conv(1:i);
end

X = reshape(X',[h w nBasis]);

end


function res = applyAt( meas, spatialSmooth, nonnegativity, h, w, cameraMat, basisFcns, beta, gamma, rho)

nBasis = size(basisFcns,2);
s1 = basisFcns'*cameraMat'*meas;

if beta > 0
    
    % Derivatives in the x direction
    spatialXSm = spatialSmooth(:,1:h*(w-1));
    spatialXSm = reshape(spatialXSm',[h, w-1, nBasis]);
    
    dx2Img = imfilter(spatialXSm,[-1 1 0]);
    dx2Img = dx2Img(:,2:w-1,:);
    dx2Img = cat(2,spatialXSm(:,1,:),dx2Img,-spatialXSm(:,w-1,:,:));
    
    dx2ImgMat = reshape(dx2Img,h*w,nBasis)';
    
    % derivatives in the y direction
    spatialYSm = spatialSmooth(:,h*(w-1)+1:h*(w-1) + w*(h-1));
    spatialYSm = reshape(spatialYSm', [h-1,w,nBasis]);
    
    dy2Img = imfilter(spatialYSm,[-1 1 0]');
    dy2Img = dy2Img(2:h-1,:,:);
    dy2Img = cat(1,spatialYSm(1,:,:),dy2Img,-spatialYSm(h-1,:,:));
    
    dy2ImgMat = reshape(dy2Img,h*w,nBasis)';
    
    s34 = rho*(dy2ImgMat + dx2ImgMat);
else
    s34 = 0;
end

if gamma > 0
    s5 = rho*basisFcns'*nonnegativity;
else
    s5 = 0;
end
res = s1 + s34 + s5;
res = res(:);

end

function res = applyAtA( meas, h,w, cameraMat, basisFcns, alpha, beta, gamma, rho)

% meas will be passed as a vector
nBasis = size(basisFcns,2);
nWaves = size(basisFcns,1);

meas = reshape(meas,nBasis,h*w);
Rlambda = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

% Measurement approximation
M1 = basisFcns'*(cameraMat'*cameraMat)*basisFcns;

% Smoothness
if alpha > 0
    M2 = alpha*basisFcns'*(Rlambda'*Rlambda)*basisFcns;
else
    M2 = 0;
end
s12 = (M1 + M2)*meas;

% Now the spatial filtering
if beta > 0
    img = reshape(meas',h,w,nBasis);
    
    % Derivatives in the x direction
    firstLastColumn = imfilter(img(:,[1 2 w-1 w],:),[1 -1 0],'same');
    interior = imfilter(img,[-1 2 -1],'same');
    interior(:,1,:) = firstLastColumn(:,2,:);
    interior(:,w,:) = -firstLastColumn(:,4,:);
    
    dx2ImgMat = reshape(interior,h*w,nBasis)';
    
    % Derivatives in the y direction
    firstLastRow = imfilter(img([1 2 h-1 h],:,:), [1 -1 0]','same');
    interior = imfilter(img,[-1 2 -1]','same');
    interior(1,:,:) = firstLastRow(2,:,:);
    interior(h,:,:) = -firstLastRow(4,:,:);
    
    dy2ImgMat = reshape(interior,h*w,nBasis)';
    
    
    s34 = rho*(dx2ImgMat + dy2ImgMat);
else
    s34 = 0;
end

% Box constraint
if gamma > 0
    M3 = rho*(basisFcns'*basisFcns);
    s5 = M3*meas;
else
    s5 = 0;
end

res = s12 + s34 + s5 ;
res = res(:);

end




