function plot_results(txIn, rxOut, config)

 [numErr, ber] = biterr(rxOut, txIn); % Calculation BER
 
% Print BER 
fprintf('For our channel BER = %f\nNumber of errors = %d\n', ...
    ber, numErr)

end