function [carrierSignal] = tx(data, config)
% Modulation
%genSinSignal = exp(2*pi*1j* config.tx.centralFrequency * config.timePoints); % generating complex sinewave  by exponent 
genSinSignal = sin(2*pi* config.tx.centralFrequency * config.timePoints); % generating sinewave NOT complex  
figure(34)
subplot(2,1,1)
plot(genSinSignal)
hold on;
grid on;
% grey coding
%data = bin2gray(data, 'psk', 2);

%Upsample
symbGrayUpsample = repelem(data, config.tx.oversampling);

txInphase   = 2*symbGrayUpsample - ones(size(symbGrayUpsample));
%txQuadrture = zeros(size(symbGrayUpsample));    % becose BPSK is not QPSK

txSignal = txInphase;   % + 1j*txQuadrture;
% figure; plot (txSignal, 'o');
% axis('equal')

carrierSignal = genSinSignal .* txSignal;

end