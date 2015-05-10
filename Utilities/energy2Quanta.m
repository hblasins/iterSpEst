function res = energy2Quanta(data,wave)

hc = 1.9865e-25;
wave = wave/1e9; %Nanometers to meters
res = data.*repmat(wave(:),[1 size(data,2)])/hc;
    
end