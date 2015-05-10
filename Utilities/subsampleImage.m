function [ imgROI ] = subsampleImage( img, h, w )

imgH = size(img,1);
imgW = size(img,2);
nChannels = size(img,3);

delta = round((imgH/h + imgW/w)/2);
deltaHalf = round(delta/2);
deltaQt = round(deltaHalf/2);

imgROI = zeros(h*deltaHalf,w*deltaHalf,nChannels);

for hh=1:h
    for ww=1:w
   
        imgROI(deltaHalf*(hh-1)+1:deltaHalf*hh,deltaHalf*(ww-1)+1:deltaHalf*ww,:) = ...
            img(deltaQt + delta*(hh-1)+1:deltaQt + delta*(hh-1) + deltaHalf,deltaQt + delta*(ww-1)+1:deltaQt + delta*(ww-1) + deltaHalf,: );
    end
end


end

