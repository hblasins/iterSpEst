function showResults( X, reflBasis, cameraMat, measVals, wave, reference )

h = size(X,1);
w = size(X,2);
nBasis = size(X,3);
nWaves = size(reflBasis,1);
nChannels = size(cameraMat,1);

reflEst = reflBasis*reshape(X,h*w,nBasis)';
reflEstImg = reshape(reflEst',[h w nWaves]);
predVals = cameraMat*reflEst;
predValsImg = reshape(predVals', [h w nChannels]);

figure; 
hold on; grid on; box on;
plot(measVals(:),predValsImg(:),'.');
xlabel('Measured');
ylabel('Predicted');
title('Camera pixel intensities');


RMSEls = sqrt(sum((reflEstImg - reference).^2,3));
figure; imagesc(RMSEls); colorbar;
title('Per pixel squared error (SE)');


dh = floor(h/4);
dw = floor(w/6);

figure;
for xx=1:6
    for yy=1:4
    
        i = 6*(yy-1) + xx;
        
        meas = reflEstImg(dh*(yy-1)+1:yy*dh,dw*(xx-1)+1:xx*dw,:);
        mr = reshape(meas,dh*dw,nWaves)';
        
        ref = reference(dh*(yy-1)+1:yy*dh,dw*(xx-1)+1:xx*dw,:);
        rr = reshape(ref,dh*dw,nWaves)';
        
        subplot(4,6,i);
        hold on; grid on; box on;
        plot(wave,mr,'g');
        plot(wave,rr,'r');
        ylim([0 1]);
        xlim([min(wave) max(wave)]);
        
    end
end



end

