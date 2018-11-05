function [rxOut] = rx(data, config)

% Modulation
%genSinSignal = exp(2*pi*1j* config.tx.centralFrequency * config.timePoints); % generating complex sinewave  by exponent 
genSinSignal = sin(2*pi* config.tx.centralFrequency * config.timePoints); % generating sinewave NOT complex  
figure(34)
subplot(2,1,2)
plot(genSinSignal, 'r')

intermedSignal = genSinSignal .* data;

%The following code is taken from MatLab Documentation
  %    Design a lowpass filter with a passband-edge frequency of 2 MHz, a 
  %    stopband-edge of 2.1 MHz, passband ripple of 0.01, stopband ripple 
  %    of 0.1, and a sampling frequency of config.tx.sampleFrequency:
 
          [n,fo,mo,w] = firpmord( [1e6 1.2e6], [1 0], [0.01 0.1], config.tx.sampleFrequency );
          b = firpm(n,fo,mo,w);
basebandSignal = filter(b, 1, intermedSignal);

% downsample
for i = 1:length(basebandSignal)/config.tx.oversampling
        
        avrBasebandSignal(i) = mean(basebandSignal((i-1)*config.tx.oversampling + 1 : i*config.tx.oversampling));
end

indexsOnes = find(avrBasebandSignal > 0.25);

onesMat = ones(size(avrBasebandSignal));
receiverSolving  = zeros(size(avrBasebandSignal));
receiverSolving(indexsOnes) = onesMat(indexsOnes);
rxOut = [receiverSolving(3:config.numBits) 0 0];

end
