function [signalWithNoise] = add_noise(data, configAwgn)

switch configAwgn.mode 
    case 1      % awgn 
        data = data + configAwgn.normDispersion * randn(size(data));
    otherwise
       
end
signalWithNoise = data.*exp(1j*configAwgn.phaseOffset); % add constant phase shift

end