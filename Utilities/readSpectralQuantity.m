function [ res ] = readSpectralQuantity( fName, wave )

value = load(fName);
res = interp1(value.wavelength,value.data,wave);


end

