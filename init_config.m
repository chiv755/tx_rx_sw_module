function config = init_config()

config.numBits = 2^8;    % Data size
config.samplesNum = 2^16;   % simulation time in smples
config.normCoeff =2;       % coefficient for normalizing amplitude

% Config for transmiter
config.tx.sampleFrequency = 160e6; % system sample rate 160MHz
config.tx.centralFrequency = 2e6;   % central frequency (mixer frequency) 2MHz
%config.tx.dataRate = 1e4;    % data rate 1Mbit/s (actually Msymbol/s)
config.tx.oversampling = config.samplesNum/config.numBits;    % oversampling rate


config.timePoints = 1/config.tx.sampleFrequency : 1/config.tx.sampleFrequency : config.samplesNum/config.tx.sampleFrequency;

% Config for channel
config.awgn.mode = 1;  % Mode additional noise
                       % 1 - awgn
                       % otherwise nothing
                       
config.awgn.normDispersion = 10;%0.5/(config.normCoeff*0.5); % normalize dispersion for amplitude noise
config.awgn.phaseOffset = 0;  % phase ofset 

% Config for reciever
%config.rx;

end